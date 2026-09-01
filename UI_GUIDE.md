# ONT Plasmid Analyzer 사용 가이드

이 UI는 외부 서버로 데이터를 전송하지 않고 **현재 Linux 워크스테이션 안에서만**
실행됩니다. 기존 `run_pipeline.sh` 명령줄 방식도 그대로 사용할 수 있습니다.

공개 저장소에는 특정 회사의 이름이나 로고가 포함되지 않습니다. 사내 CI/브랜드
담당자가 승인한 설정은 프로젝트 루트의 `branding.local.json`과
`branding_logo.svg`(또는 PNG)에 넣을 수 있으며, 이 파일들은 Git에서 자동 제외됩니다.
왼쪽 하단 로고는 `branding_sidebar_logo.svg` 또는 `branding_sidebar_logo.png`로 넣습니다.
`branding.example.json`을 복사해 조직명, 문구 및 색상을 설정할 수 있습니다.

## 1. 설치

기존 `NGS_env`를 업데이트합니다.

```bash
cd NGS_ONT
conda env update -p /home/MCET03/conda_envs/NGS_env -f environment.yml
conda activate /home/MCET03/conda_envs/NGS_env
```

인터넷이 차단된 사내 서버에서는 로컬 미러만 지정할 수 있습니다.

```bash
conda env update -p /home/MCET03/conda_envs/NGS_env --offline --override-channels \
  -c file:///data/conda_repo/mirror/conda-forge/conda-forge \
  -c file:///data/conda_repo/mirror/bioconda \
  -f environment.yml
```

설치 확인:

```bash
streamlit --version
minimap2 --version
samtools --version
bcftools --version
NanoFilt --version
```

## 2. UI 실행

```bash
conda activate /home/MCET03/conda_envs/NGS_env
cd /data/user/MCET03/04_ONT/NGS_ONT
./run_ui.sh
```

브라우저에서 다음 주소를 엽니다.

```text
http://<Linux-server-IP>:8501
```

포트를 바꾸려면:

```bash
ONT_UI_PORT=8502 ./run_ui.sh
```

서버 외부 접속을 막고 SSH 포트 포워딩만 사용하려면:

```bash
ONT_UI_ADDRESS=127.0.0.1 ./run_ui.sh
```

## 3. Batch analysis (권장)

1. 왼쪽 `분석 메뉴`에서 **Batch analysis**를 선택합니다.
2. 필요한 reference FASTA를 첫 번째 업로드 영역에 모두 드래그합니다.
3. Reference는 파일명 기준 자연 정렬되며(예: `01, 02, ..., 10`),
   `Reference file check`에서 정렬된 파일명과 길이를 확인합니다.
4. Reference 파일명별로 생성된 영역에 해당 FASTQ 파일 또는 sample 폴더를 그대로
   드래그합니다. 별도의 입력 형식 선택은 필요하지 않습니다.
5. 필요하면 `Sample ID 확인/수정`을 열어 자동 인식된 ONT sample name을 수정합니다.
   여러 FASTQ chunk에 같은 Sample ID를 입력하면 하나의 sample로 합쳐 분석합니다.
6. `분석 전 최종 확인`에서 reference별 Sample ID와 FASTQ 개수를 확인합니다.
7. 왼쪽의 `Analysis settings` 버튼을 누르면 중앙에 나타나는 `Batch settings`에서
   CPU, 병렬 sample 수, read filter 및 variant review 기준을 저장할 수 있습니다.
8. 중복·누락 경고가 없는지 확인한 뒤 `Batch analysis 실행`을 누릅니다.

각 입력 항목의 `?`에 커서를 올리면 설정 의미와 권장 사용법을 확인할 수 있습니다.

결과는 reference별로 연결된 sample이 묶여 표시되며 `CLEAN`,
`VARIANT DETECTED`, `REVIEW`, `ERROR` 상태를 제공합니다. 전체 결과 요약은 CSV로
내려받을 수 있고 BAM, VCF, consensus FASTA 및 보고서는 서버 결과 폴더에 남습니다.
실행 command와 상세 출력은 기본적으로 접힌 `Analysis log`에서 필요할 때만 확인합니다.

### 내장 예제 데이터

Batch 화면의 `예제 데이터로 테스트하기`를 열어 ZIP을 받습니다.
압축을 푼 뒤 `references`의 FASTA 5개를 올립니다. 각 reference 입력창에
`demo_reads/<같은 reference 이름>/`의 FASTQ 파일 또는 폴더를 그대로 올립니다.
예제는 reference마다 exact/SNP/insertion sample을
하나씩 포함하며 `expected_mapping.csv`에서 예상 연결과 변이를 확인할 수 있습니다.

Sample ID는 `barcode13` 같은 번호, ONT sample alias 폴더명 또는 FASTQ 파일명에서
자동으로 가져옵니다. 따라서 barcode라는 이름을 반드시 사용할 필요가 없습니다.

> 디렉터리 업로드를 위해 Streamlit 1.57 이상이 필요합니다. FASTQ는 브라우저를
> 통해 서버로 전송되므로 분석 중 브라우저 탭을 닫지 마세요.

## 4. Quick FASTA comparison

왼쪽 `분석 메뉴`에서 **Quick sequence comparison**을 선택한 뒤 reference와
query/consensus를 각각 FASTA 업로드 또는 서열 붙여넣기로 입력합니다.

- **Reference**: 원래 설계한 plasmid/vector의 기준 FASTA입니다.
- **Query**: ONT 분석 후 얻은 consensus FASTA, assembly 결과 또는 확인하려는 완성
  plasmid 서열입니다. Raw FASTQ는 Query에 넣지 않습니다.
- Raw ONT read부터 분석하려면 **Batch analysis**를 사용합니다.

- 방향과 reverse complement를 자동 판별합니다.
- `Circular plasmid/vector`를 켜면 FASTA 시작/끝 이음매를 통과하는 정렬을 허용합니다.
- SNP(1-bp point mutation), insertion, deletion의 위치와 주변 reference context를 표시합니다.
- 결과 표는 CSV로 받을 수 있습니다.

이 모드는 **두 완성 서열 사이의 차이**만 보여주므로 read depth, QUAL, allele
fraction은 제공하지 않습니다. 해당 신뢰도 지표가 필요하면 Raw ONT analysis를
사용합니다.

## 5. 모드 선택 기준

- **Batch analysis**: raw ONT FASTQ를 분석할 때 사용합니다. Reference 하나와 해당
  sample folder 하나만 올리면 single-sample 분석도 같은 화면에서 실행할 수 있습니다.
- **Quick sequence comparison**: 이미 완성된 query/consensus FASTA가 있고
  reference와 sequence 차이만 빠르게 확인할 때 사용합니다. Read depth, QUAL,
  allele fraction은 계산하지 않습니다.
- 기존 **Single-sample ONT analysis** 코드는 보존하지만 일반 메뉴에서는 숨깁니다.
  Batch analysis와 기능이 겹쳐 부서원용 화면을 불필요하게 복잡하게 만들기 때문입니다.

## 6. 결과 위치

각 실행은 서로 섞이지 않게 다음 경로에 저장됩니다.

```text
ui_runs/<실행시각>_<experiment>_<job-id>/
├── config.yaml
├── pipeline.log
├── references/
├── data/          # 서버 FASTQ는 복사하지 않고 symbolic link로 연결
└── results/
```

UI에서 CSV, VCF.GZ, consensus FASTA, sample report를 내려받을 수 있습니다.
BAM/BAM.BAI와 depth 파일은 크기가 클 수 있어 결과 폴더에만 보관합니다.

## 7. 정확도 해석

- SNP는 이 UI에서 1-bp substitution/point mutation을 의미합니다.
- `PASS`는 현재 선택한 QUAL, DP, allele-fraction, edge 기준을 만족합니다.
- `REVIEW`는 오류 확정이 아니라 추가 확인이 필요한 call입니다.
- ONT indel은 homopolymer에서 false positive가 증가할 수 있으므로 BAM을 IGV로
  열어 read-level support를 함께 확인하는 것이 좋습니다.
- plasmid/vector 단일 clone의 진짜 변이는 보통 높은 allele fraction을 보입니다.

## 8. 테스트

```bash
python -m unittest discover -s tests -v
bash -n run_pipeline.sh run_ui.sh
```
