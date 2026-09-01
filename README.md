# Nanopore (MinION) Reference-Mapping Pipeline

Oxford Nanopore MinION으로 시퀀싱한 결과(fastq)를 reference FASTA(벡터맵)에
매핑하고, BAM/consensus FASTA/VCF/리포트를 생성하는 파이프라인입니다.

## 로컬 웹 UI: ONT Plasmid Analyzer

명령어 없이 브라우저에서 분석할 수 있는 오프라인 Streamlit UI가 포함되어
있습니다. 데이터와 분석은 Linux 서버 내부에서만 처리되며 Windows 사용자는
Chrome/Edge로 서버 주소에 접속합니다.

```bash
conda env update -n NGS_env -f environment.yml
conda activate NGS_env
./run_ui.sh
```

브라우저에서 `http://<Linux-server-IP>:8501`을 열면 왼쪽 `분석 메뉴`에서
다음 두 모드를 선택할 수 있습니다. `Analysis settings` 버튼을 누르면 선택한
분석의 설정이 중앙 작업공간에 표시되며, 각 `?` 도움말에 커서를 올리면 설명을
확인할 수 있습니다. 앱은 일관된 light theme와 브랜드 그라데이션 sidebar를 사용합니다.

1. **Batch analysis**: 이번 run의 전체 barcode sample 수를 1–96개 중
   선택하고 필요한 reference와 ONT barcode 폴더를 업로드합니다. 각 barcode를
   원하는 reference에 독립적으로 연결하므로 reference별 sample 수는 자유롭습니다.
   화면에서 바로 받을 수 있는 5-reference/15-sample synthetic demo도 제공합니다.
2. **Quick sequence comparison**: reference와 query/consensus 서열을 붙여넣거나
   업로드하여 SNP(point mutation), insertion, deletion을 즉시 비교합니다.
   reverse complement와 circular plasmid 이음매를 자동으로 처리합니다.

기존 **Single-sample ONT analysis** 기능은 코드에 보존하지만 일반 메뉴에서는
숨깁니다. Batch analysis를 sample 수 1개로 설정하면 동일한 일상 분석을 더 단순한
화면에서 수행할 수 있기 때문입니다.

대용량 FASTQ는 브라우저로 업로드하지 않고 서버 경로를 선택하는 방식을
권장합니다. 기본 read QC는 `500 bp / Q10`이며 coverage가 부족할 때만
`300 bp / Q8` 완화를 고려합니다. 설치·화면·결과 해석은
[`UI_GUIDE.md`](UI_GUIDE.md)를 참고하세요.


## 핵심 아이디어: 폴더명 기반 자동 매칭

`references/`와 `data/` 아래에 **"날짜_실험명" 폴더를 같은 이름으로** 만들어
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
├── data/
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

- `references/<실험폴더>/` 안에 reference fasta가 **여러 개** 있으면, 각
  reference마다 매칭되는 fastq를 찾아 결과 폴더를 따로 만듭니다.
- `references/`에는 있지만 `data/`에 같은 이름 폴더가 없으면(또는 반대)
  해당 실험은 건너뛰고 경고만 출력합니다.

### reference ↔ fastq 매칭 우선순위

각 reference fasta(`<reference_name>.fasta`)에 대해, 다음 순서로 매칭되는
fastq를 찾습니다:

1. **파일명이 같은 파일**: `data/<실험폴더>/` 아래(하위 폴더 포함) 어디에든
   `<reference_name>.*` (예: `R34.141-DAR-revSPG23_TIR90.gz`,
   `R34.141-DAR-revSPG23_TIR90.fastq.gz` 등 확장자 무관)인 파일이 있으면
   그 파일을 사용합니다. reference 파일명과 데이터 파일명(확장자 제외)이
   완전히 동일한 경우에 적용됩니다.
2. **이름이 같은 폴더**: `data/<실험폴더>/<reference_name>/` 가 있으면
   그 안의 fastq를 사용 (sample sheet의 alias로 demultiplex한 경우).
3. **번호 매칭**: reference 파일명이 `01_151-NIV-fwdGS_TIR68.fasta`처럼
   **앞자리 숫자**로 시작하면, 그 숫자와 같은 번호의
   `data/<실험폴더>/barcode01/`(barcode + 같은 숫자) 폴더를 자동으로
   찾아 사용합니다. 즉 96-well에 reference를 `01_..., 02_..., ... 96_...`
   처럼 번호를 붙여두면, 시퀀싱 결과의 `barcode01`~`barcode96`과 자동으로
   1:1 매핑됩니다 (앞자리 0 유무는 무시: `01`과 `1` 모두 `barcode01`과
   매칭).
4. **공유 폴백**: 위 세 가지로 못 찾으면, 실험 폴더 전체의 fastq(모든
   reference/barcode 서브폴더 및 위에서 매칭된 파일 제외)를 사용합니다 —
   reference가 1개뿐인 실험에 적합합니다.

> 데이터 파일은 `.fastq`, `.fastq.gz`, `.fq`, `.fq.gz` 뿐 아니라 확장자가
> `.gz`만 있는 경우(예: `R34.141-DAR-revSPG23_TIR90.gz`)도 인식합니다.
> reference 파일에 `.fa.amb`, `.fa.ann`, `.fa.bwt`, `.fa.pac`, `.fa.sa`
> 같은 BWA 인덱스 파일이 같이 있어도 무시되며 (`.fasta`/`.fa`/`.fna`만
> reference로 사용), 분석에 영향을 주지 않습니다.

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
별도 작업 없이 바로 `data/<실험폴더>/`에 옮기고 `run_pipeline.sh`만
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
5. 그 출력 폴더 전체를 `data/<날짜>_<실험명>/` (references와 같은
   폴더명) 아래로 옮기고 `./run_pipeline.sh`를 실행한다.

> 이미 barcode 번호로 시퀀싱이 끝난 데이터가 있고, 그 실험 폴더 안에
> reference가 **1개뿐**이라면 이 단계는 건너뛰어도 됩니다 —
> `data/<날짜>_<실험명>/barcode01/...` 형태 그대로 두면
> `run_pipeline.sh`가 그 fastq를 해당 reference 1개에 매핑합니다.
> reference가 **여러 개**인 실험이라면, 각 reference의 fastq가 어느
> 것인지 구분이 필요하므로 위 sample sheet 방식(또는 수동으로
> `data/<실험폴더>/<reference이름>/`에 해당 fastq를 넣는 방식)을
> 사용해야 합니다.

## 사전 준비

1. 필요한 도구 설치 (conda 권장):
   ```bash
   conda install -c bioconda -c conda-forge minimap2 samtools bcftools nanofilt
   ```
2. `references/<날짜>_<실험명>/`에 벡터맵 fasta 파일을 넣는다.
3. `data/<날짜>_<실험명>/`에 (위와 동일한 폴더명으로) MinKNOW/Guppy/Dorado의
   barcode별(또는 alias별) fastq.gz 폴더를 넣는다.
4. 필요하면 `config.yaml`에서 minimap2 preset, threads, QC 필터링 기준,
   variant caller(`bcftools`, `medaka`, 또는 `pilon`)를 조정한다.
   `pilon`을 쓰려면 별도로 설치된 Java와 `pilon-*.jar` 경로
   (`pilon_jar`)가 필요하다 (자세한 내용은 아래 "Pilon 기반 변이/오류
   검출" 참고).

## 실행

```bash
cd nanopore_pipeline
./run_pipeline.sh
```

`references/`와 `data/` 아래에서 **이름이 일치하는 폴더 쌍**을 모두
찾아서 번호 목록으로 보여줍니다. 새 실험을 추가할 때는 두 폴더 아래에
같은 이름의 폴더만 만들어주면 자동으로 목록에 나타납니다.

```
Available experiments:
   1) 260609_newexp
   2) 260610_newexp

Run which experiment(s)? (number(s), comma-separated, or Enter for all):
```

- 번호를 입력하면 해당 실험만 실행 (예: `1` 또는 `1,2`)
- 그냥 Enter를 누르면 전체 실험 실행
- 터미널이 아닌 환경(cron 등)에서 실행하면 입력 없이 자동으로 전체 실행

### 결과 폴더 이름 (실행 시점 타임스탬프 포함)

결과는 `results/<실험명>_<실행한 날짜>_<시간>/<reference_name>/`에
저장됩니다. 예를 들어 `260609_exp` 실험을 `2026-06-10 18:10`에
실행하면:

```
results/260609_exp_260610_1810/<reference_name>/
```

같은 실험을 다시 실행하면 새로운 타임스탬프 폴더가 생성되어, 이전
실행 결과와 섞이지 않고 따로 보관됩니다.

## 병렬 처리 (고성능 워크스테이션)

각 reference(샘플)는 서로 독립적으로 처리되므로, `config.yaml`의
`parallel_jobs`를 올리면 여러 샘플을 동시에 처리할 수 있습니다.

- `threads`: 샘플 1개당 minimap2/samtools가 사용할 CPU 코어 수
- `parallel_jobs`: 동시에 처리할 샘플(reference) 개수

`threads * parallel_jobs`가 워크스테이션의 CPU 코어 수를 넘지 않도록
설정하세요. 예를 들어 **256코어** 워크스테이션이라면:

```yaml
threads: 8
parallel_jobs: 32
```

처럼 설정하면 8코어 × 32개 샘플 = 256코어를 모두 활용해서 모든
샘플을 동시에 처리합니다. 샘플 수가 32개보다 적으면 그만큼만 병렬로
돌고, 더 많으면 먼저 끝난 샘플부터 다음 샘플을 이어서 처리합니다.

## 전체 결과 요약 (변이 한눈에 보기)

샘플이 많을 때 `*_report.md`/`*.vcf.gz`를 하나씩 열어보지 않고, 모든
실험/샘플의 매핑률·평균 depth·변이(point mutation/insertion/deletion)를
표 하나로 정리할 수 있습니다:

```bash
python scripts/summarize_variants.py
```

실행하면 `results/` 아래의 결과 폴더 목록을 보여주고, 어떤 결과를
요약할지 물어봅니다:

```
Available results:
   1) 260609_exp_260610_1810
   2) 260610_newexp_260610_1900

Summarize which result(s)? (number(s), comma-separated, or Enter for all):
```

- 번호 하나를 선택하면 (예: `1`), 그 결과 폴더와 같은 이름의 CSV가
  `results/260609_exp_260610_1810.csv`로 생성됩니다.
- 여러 개(`1,2`) 또는 Enter(전체)를 선택하면 `results/summary_report.csv`로
  생성됩니다.
- CSV는 Excel에서 바로 열어서 정렬/필터할 수 있습니다. Markdown 표로
  받고 싶으면 `--output results/summary_report.md`처럼 확장자를 `.md`로
  지정하세요.


### 노이즈성 변이 필터링 (`--min-qual`, `--min-depth`)

ONT 데이터는 에러율이 높아서, 특히 homopolymer(같은 염기 반복) 구간에서
`bcftools` 기본 설정으로는 실제로는 시퀀싱 에러인 indel/SNP까지 변이로
잡히는 경우가 많습니다. `QUAL`(변이 신뢰도)과 `DP`(depth)에 최소 기준을
줘서 신뢰도 낮은 변이를 제외할 수 있습니다:

```bash
python scripts/summarize_variants.py --min-qual 20 --min-depth 10
```

- `--min-qual 20`: VCF의 QUAL 컬럼이 20 미만인 변이는 제외
- `--min-depth 10`: 해당 위치의 read depth(`INFO/DP`)가 10 미만인 변이는 제외

기준값은 데이터 특성에 맞게 조정하세요. 필터링 후에도 의심되는 변이는
`bcftools view`로 원본 VCF를 직접 열어 QUAL/DP를 확인하거나, IGV로
`.sorted.bam`을 열어 해당 위치의 read를 직접 봐서 진짜 변이인지
판단하는 것이 좋습니다.

### 레퍼런스 양 끝(이음매) 변이 분리 (`--edge-margin`)

원형 플라스미드/벡터를 선형 fasta로 펼치면, 실제로는 이어져 있는
부위가 레퍼런스 상에서는 양쪽 끝(1번 위치 근처, 마지막 위치 근처)으로
나뉘게 됩니다. 이 "이음매"를 가로지르는 read는 정렬이 어긋나기 쉬워서,
레퍼런스 양 끝 부근에 같은 패턴의 indel이 반복적으로 잘못 호출되는
경우가 많습니다.

결과 선택 후, 다음과 같이 추가로 물어봅니다:

```
Separate variants near reference ends (circular-plasmid junction artifacts)? Edge margin in bp [default 50, 0 = off]:
```

- 그냥 Enter → 기본값 50bp 적용
- `0` → 끝부분 분리 안 함 (기존 동작과 동일)
- 원하는 숫자 입력 → 레퍼런스 시작/끝 그만큼(bp) 안에 있는 변이를
  `Variants`/`SNP`/`Ins`/`Del` 집계에서 제외하고 별도의
  `Edge variants` 컬럼에 표시

ONT 데이터에서는 보통 **50bp** 정도면 충분하지만, depth가 매우 높거나
read가 짧다면 100bp까지 늘려도 됩니다.

`Edge variants`에 표시된 항목은 진짜 변이가 아닌 매핑 아티팩트일
가능성이 높지만, 혹시 실제 벡터 구조상 그 위치에 진짜 차이가 있을
수도 있으니 한 번씩 확인하는 것을 권장합니다.

명령줄에서 바로 지정하고 싶으면 (non-interactive) `--edge-margin 50`
옵션을 줄 수도 있습니다.

### ONT 전용 변이 호출 설정 (`-X ont`, `--ploidy 1`)

`run_pipeline.sh`는 `bcftools mpileup`이 `-X`/`--config` 프리셋을
지원하면 자동으로 `-X ont`를 적용합니다. 이 옵션은 ONT 리드의 indel
에러 프로파일(특히 homopolymer 구간)을 고려해서 indel/SNP 호출
파라미터를 조정해주므로, 기본값보다 false-positive 변이가 줄어듭니다.

또한 `bcftools call`에는 `--ploidy 1`을 지정합니다. 플라스미드/벡터는
단일 클론(haploid)이므로, 진짜 변이라면 거의 모든 read(>90%)가 그 변이를
지지해야 합니다. 기본값(diploid)으로는 read의 약 절반 정도만 지지하는
homopolymer 에러도 `0/1` heterozygous 변이로 호출되어 노이즈가 많이
잡히는데, `--ploidy 1`을 쓰면 0/0(레퍼런스) 또는 1/1(변이) 중 하나만
선택하게 되어 이런 애매한 노이즈가 변이로 호출될 가능성이 크게
줄어듭니다.

기존 결과에 다시 적용하려면 그냥 `./run_pipeline.sh`를 다시 실행하면
됩니다 (read 병합 단계는 결과가 이미 있으면 건너뛰고, mapping/consensus/
variant calling은 새 옵션으로 다시 생성됩니다). 그 후
`summarize_variants.py`도 다시 실행해서 요약을 갱신하세요.

### Pilon 기반 변이/오류 검출 (`variant_caller: pilon`)

기존에 사용하던 BWA + Pilon 기반 스크립트와 동일한 방식으로 변이를
검출하고 싶다면 `config.yaml`에서 다음과 같이 설정하세요:

```yaml
variant_caller: pilon
pilon_jar: /home/mcet/install_file/usbmnt/17.pilon/pilon-1.24.jar
pilon_mem: 16G
```

- `pilon_jar`: `pilon-*.jar`의 절대 경로 (또는 `config.yaml`이 있는
  폴더 기준 상대 경로). 이 값이 비어있거나 파일을 찾을 수 없으면 해당
  샘플의 variant calling 단계는 건너뜁니다 (경고만 출력).
- `pilon_mem`: Pilon(Java)에 할당할 최대 힙 메모리 (`-Xmx` 값). 워크스테이션의
  가용 메모리에 맞게 조정하세요 (예: `200G`).
- Pilon은 conda 환경이 아닌 별도로 설치된 Java(JRE/JDK)가 필요합니다.
  `java -version`이 동작하는지 먼저 확인하세요.

`variant_caller: pilon`일 때 5단계(variant calling)는 minimap2로 만든
정렬 결과(`<reference_name>.sorted.bam`)에 대해 다음 명령으로 실행됩니다:

```bash
java -Xmx<pilon_mem> -jar <pilon_jar> \
    --genome <reference>.fasta --bam <reference_name>.sorted.bam \
    --output <reference_name> --outdir results/.../pilon \
    --fix all --mindepth 2.0 --changes --vcf --verbose --threads <threads>
```

생성되는 파일:

- `<reference_name>.changes`: Pilon이 실제로 수정(교정)한 위치 목록
  (`old_locus new_locus old_seq new_seq` 형식). 결과 폴더 최상위에도
  복사되어 `summarize_variants.py`가 자동으로 읽습니다.
- `<reference_name>.vcf.gz`: Pilon의 `--vcf` 출력 중 변화가 없는 위치
  (`ALT="."`)는 제거하고, 실제 변이/교정 레코드만 남긴 VCF
  (다른 variant caller와 동일한 형식으로 맞춰서 저장).
- 원본 Pilon 출력(`*.changes`, `*.vcf`, `*.fasta`, 로그 등)은
  `results/.../<reference_name>/pilon/` 폴더에 그대로 남아있습니다.
