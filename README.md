# Nanopore (MinION) Reference-Mapping Pipeline

Oxford Nanopore MinION으로 시퀀싱한 결과(fastq)를 reference FASTA(벡터맵)에
매핑하고, BAM/consensus FASTA/VCF/리포트를 생성하는 파이프라인입니다.

## 핵심 아이디어: 폴더명 기반 자동 매칭

`references/`와 `data/raw/` 아래에 **"날짜_실험명" 폴더를 같은 이름으로** 만들어
두면, 파이프라인이 그 둘을 자동으로 매칭해서 매핑 분석을 진행합니다.
출력 파일은 항상 **reference FASTA 파일명(실제 벡터/샘플 이름)** 기준으로
생성되므로, barcode 번호 같은 시퀀서 이름을 분석 후에 따로 바꿔줄 필요가
없습니다.

```
nanopore_pipeline/
├── config.yaml
├── run_pipeline.sh
├── references/
│   └── 20260610_pUC19_test/        <- 실험 폴더 (날짜_실험명)
│       └── pUC19_insertA.fasta     <- 벡터맵 전체 서열 (이게 "real name")
├── data/raw/
│   └── 20260610_pUC19_test/        <- 같은 이름의 폴더
│       └── barcode01/              <- MinKNOW/Guppy/Dorado 결과 (서브폴더 있어도 됨)
│           └── *.fastq.gz
└── results/
    └── 20260610_pUC19_test/
        └── pUC19_insertA/           <- barcode01이 아니라 reference 이름으로 생성됨
            ├── pUC19_insertA.sorted.bam / .bai
            ├── pUC19_insertA.flagstat.txt
            ├── pUC19_insertA.depth.txt
            ├── pUC19_insertA.consensus.fasta
            ├── pUC19_insertA.vcf.gz
            └── pUC19_insertA_report.md
```

- `references/<실험폴더>/` 안에 reference fasta가 **여러 개** 있으면, 그
  실험 폴더의 fastq 전체를 각 reference에 대해 모두 매핑해서 reference별로
  결과 폴더를 따로 만듭니다.
- `data/raw/<실험폴더>/` 안의 fastq.gz는 barcode 서브폴더에 있든 바로 있든
  상관없이 재귀적으로 모두 찾아서 합칩니다.
- `references/`에는 있지만 `data/raw/`에 같은 이름 폴더가 없으면(또는 반대)
  해당 실험은 건너뛰고 경고만 출력합니다.

## 레퍼런스 FASTA 관련 Q&A

**Q. 레퍼런스에는 읽고 싶은 벡터맵 전체를 넣으면 되나요?**
네. 플라스미드/벡터의 전체 서열 1개를 FASTA 파일 하나에 넣으면 됩니다
(원형이면 임의의 한 지점에서 잘라 선형 서열로 저장).

**Q. Forward/Reverse를 나눠서 레퍼런스를 따로 만들어야 하나요?**
아니요. minimap2는 reference의 양쪽 가닥을 모두 검사해서 정렬하므로
정방향 서열 1개만 있으면 충분합니다.

**Q. `.fa.amb`, `.fa.ann`, `.fa.bwt`, `.fa.pac`, `.fa.sa` 같은 인덱스 파일이
필요한가요?**
아니요. 그 파일들은 **BWA**용 인덱스입니다. 이 파이프라인은 ONT 데이터에
적합한 **minimap2**를 사용하며, minimap2는 `.fasta` 원본 파일만으로 즉시
인덱싱(메모리 내) 후 매핑하므로 별도 인덱스 파일을 만들 필요가 없습니다.
(반복 실행 속도를 높이고 싶다면 `minimap2 -d ref.mmi ref.fasta`로 `.mmi`
인덱스를 미리 만들어 둘 수는 있지만, 필수는 아닙니다.)

## 시퀀싱 시작 "전" 단계: 96-well sample sheet 자동 생성 (권장)

핵심 요약: **시퀀싱 시작 전에, 각 well(barcode)에 reference 파일명과 동일한
alias를 매핑한 sample sheet를 MinKNOW에 등록해두면, basecalling/demultiplexing
결과 폴더가 처음부터 reference 이름으로 생성됩니다.** 그러면 이름을 맞추는
별도 작업 없이 바로 `data/raw/<실험폴더>/`에 옮기고 `run_pipeline.sh`만
실행하면 됩니다.

1. `references/<날짜>_<실험명>/`에 96-well 각각에 해당하는 reference fasta를
   넣는다 (파일명 = 원하는 최종 샘플명).
2. sample sheet를 자동 생성한다:
   ```bash
   python scripts/generate_samplesheet.py \
       --references references/20260610_pUC19_test \
       --experiment-id 20260610_pUC19_test \
       --output references/20260610_pUC19_test/sample_sheet.csv
   ```
   - 기본적으로 `references/<실험폴더>/` 안의 fasta 파일을 **이름순**으로
     barcode01, barcode02, ... 에 순서대로 배정합니다.
   - plate 배치 순서가 다르면 `--order order.txt` (한 줄에 reference 이름
     하나씩, plate 순서대로)로 직접 지정할 수 있습니다.
3. 생성된 `sample_sheet.csv`를 MinKNOW 실행 화면의 **Start run -> Sample
   sheet -> Browse**에서 불러온 뒤 시퀀싱을 시작한다.
4. Basecalling/demultiplexing이 끝나면 출력 폴더가 `barcode01` 대신
   sample sheet의 `alias`(= reference 파일명)로 생성된다.
5. 그 출력 폴더 전체를 `data/raw/<날짜>_<실험명>/` (references와 같은
   폴더명) 아래로 옮기고 `./run_pipeline.sh`를 실행한다.

> 이미 barcode 번호로 시퀀싱이 끝난 데이터가 있고, 그 실험 폴더 안에
> reference가 **1개뿐**이라면 이 단계는 건너뛰어도 됩니다 —
> `data/raw/<날짜>_<실험명>/barcode01/...` 형태 그대로 두면
> `run_pipeline.sh`가 그 fastq를 해당 reference 1개에 매핑합니다.
> reference가 **여러 개**인 실험이라면, 각 reference의 fastq가 어느
> 것인지 구분이 필요하므로 위 sample sheet 방식(또는 수동으로
> `data/raw/<실험폴더>/<reference이름>/`에 해당 fastq를 넣는 방식)을
> 사용해야 합니다.

## 사전 준비

1. 필요한 도구 설치 (conda 권장):
   ```bash
   conda install -c bioconda -c conda-forge minimap2 samtools bcftools nanofilt
   ```
2. `references/<날짜>_<실험명>/`에 벡터맵 fasta 파일을 넣는다.
3. `data/raw/<날짜>_<실험명>/`에 (위와 동일한 폴더명으로) MinKNOW/Guppy/Dorado의
   barcode별(또는 alias별) fastq.gz 폴더를 넣는다.
4. 필요하면 `config.yaml`에서 minimap2 preset, threads, QC 필터링 기준,
   variant caller(`bcftools` 또는 `medaka`)를 조정한다.

## 실행

```bash
cd nanopore_pipeline
./run_pipeline.sh
```

`references/`와 `data/raw/` 아래의 모든 실험 폴더를 자동으로 스캔해서,
이름이 일치하는 쌍에 대해서만 매핑을 수행합니다. 새 실험을 추가할 때는
두 폴더 아래에 같은 이름의 폴더만 만들어주면 됩니다.
