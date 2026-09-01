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

## 3. Batch plasmid analysis (권장)

1. 왼쪽 `Analysis menu`에서 **Batch plasmid analysis**를 선택합니다.
2. `Number of references`에서 이번에 분석할 plasmid 수(1–32개)를 지정합니다.
3. 지정한 수만큼 reference FASTA를 첫 번째 업로드 영역에 드래그합니다.
4. `barcode13`, `barcode14` 같은 폴더들이 들어 있는 ONT run 상위 폴더를 두 번째
   업로드 영역에서 선택합니다.
5. UI는 barcode 번호를 숫자순으로 정렬하고 reference당 3개씩 자동 배정합니다.
6. 매핑 표의 `Order` 또는 `Barcode 1–3`을 수정해 실제 실험 순서와 맞춥니다.
7. 왼쪽 `Analysis settings`에서 CPU, 병렬 sample 수, read filter를 확인합니다.
8. 중복·누락 경고가 없는지 확인한 뒤 `Run batch analysis`를 누릅니다.

각 입력 항목의 `?`에 커서를 올리면 설정 의미와 권장 사용법을 확인할 수 있습니다.

결과는 reference별로 세 replicate가 한 묶음으로 표시되며 `CLEAN`,
`VARIANT DETECTED`, `REVIEW`, `ERROR` 상태를 제공합니다. 전체 결과 요약은 CSV로
내려받을 수 있고 BAM, VCF, consensus FASTA 및 보고서는 서버 결과 폴더에 남습니다.

> 디렉터리 업로드를 위해 Streamlit 1.57 이상이 필요합니다. FASTQ는 브라우저를
> 통해 서버로 전송되므로 분석 중 브라우저 탭을 닫지 마세요.

## 4. Quick FASTA comparison

왼쪽 `Analysis menu`에서 **Quick sequence comparison**을 선택한 뒤 reference와
query/consensus를 각각 FASTA 업로드 또는 서열 붙여넣기로 입력합니다.

- 방향과 reverse complement를 자동 판별합니다.
- `Circular plasmid/vector`를 켜면 FASTA 시작/끝 이음매를 통과하는 정렬을 허용합니다.
- SNP(1-bp point mutation), insertion, deletion의 위치와 주변 reference context를 표시합니다.
- 결과 표는 CSV로 받을 수 있습니다.

이 모드는 **두 완성 서열 사이의 차이**만 보여주므로 read depth, QUAL, allele
fraction은 제공하지 않습니다. 해당 신뢰도 지표가 필요하면 Raw ONT analysis를
사용합니다.

## 5. Single-sample ONT analysis

1. 왼쪽 `Analysis menu`에서 **Single-sample ONT analysis**를 선택합니다.
2. 왼쪽 `Analysis settings`에서 read filter, caller 및 variant 기준을 설정합니다.
3. Experiment 이름과 sample/vector 이름을 입력합니다.
4. Reference FASTA를 올리거나 전체 서열을 붙여넣습니다.
5. 대용량 데이터는 `Server folder or file path`를 선택하고
   `/data/.../barcode01` 같은 경로를 입력합니다.
6. 작은 파일만 브라우저 업로드를 사용합니다.
7. `Run ONT analysis`를 누릅니다.

기본 QC는 다음과 같습니다.

- minimum read length: **500 bp**
- minimum mean read quality: **Q10**
- coverage가 부족할 때 완화 후보: **300 bp / Q8**

기본 variant review 기준은 `QUAL ≥20`, `DP ≥10`, alternate allele fraction
`≥0.80`입니다. 기준 미달 call을 삭제하지 않고 `REVIEW`로 표시합니다.
원형 reference 양 끝 50 bp와 homopolymer 인접 call도 경고로 표시됩니다.

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
