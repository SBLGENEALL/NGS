# ONT Plasmid Analyzer 사용 가이드

이 UI는 외부 서버로 데이터를 전송하지 않고 **현재 Linux 워크스테이션 안에서만**
실행됩니다. 기존 `run_pipeline.sh` 명령줄 방식도 그대로 사용할 수 있습니다.

공개 저장소에는 특정 회사의 이미지 로고를 포함하지 않습니다. 왼쪽 하단에 배경 없는
text wordmark를 사용하려면 Linux server에서 다음을 한 번 실행합니다. 이 설정 파일은
Git에서 자동 제외되며 NASCA로 암호화된 이미지가 필요하지 않습니다.

```bash
./configure_local_branding.sh "SAMSUNG BIOLOGICS" "Jongin Baek"
```

## 1. 설치

기존 `NGS_ONT_env`를 업데이트합니다.

```bash
cd NGS_ONT
conda env update -p /home/MCET03/conda_envs/NGS_ONT_env -f environment.yml
conda activate /home/MCET03/conda_envs/NGS_ONT_env
```

인터넷이 차단된 사내 서버에서는 로컬 미러만 지정할 수 있습니다.

```bash
conda env update -p /home/MCET03/conda_envs/NGS_ONT_env --offline --override-channels \
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

### 아이콘으로 실행

Linux desktop에서는 최초 한 번만 다음을 실행합니다.

```bash
cd /data/user/MCET03/04_ONT/NGS_ONT
chmod +x launch_ui.sh install_linux_launcher.sh run_ui.sh
./install_linux_launcher.sh
```

이후 바탕화면의 **ONT Plasmid Analyzer** 아이콘을 클릭하면 `NGS_ONT_env` 탐색,
UI 실행 및 브라우저 열기가 자동으로 진행됩니다. 기본 포트는 `8502`입니다.

Windows + WSL 환경에서는 프로젝트의 `Launch_ONT_UI.cmd`를 바탕화면으로 복사한 뒤
더블클릭합니다. 기본 WSL 프로젝트 경로는 `~/NGS_ONT_batch`입니다.

Windows에서 회사 Linux server에 SSH로 접속하는 환경은 server에서 다음처럼 실행합니다.

```bash
chmod +x create_remote_launcher.sh
./create_remote_launcher.sh <Linux-server-IP>
```

생성된 `Launch_ONT_UI_Remote.cmd`를 Windows 바탕화면으로 내려받아 더블클릭하면
SSH tunnel, `NGS_ONT_env`, UI 및 `localhost:8502` 브라우저가 순서대로 실행됩니다.
열린 SSH 창은 분석 중 닫지 않아야 합니다.

### 명령어로 실행

```bash
conda activate /home/MCET03/conda_envs/NGS_ONT_env
cd /data/user/MCET03/04_ONT/NGS_ONT
./run_ui.sh
```

코드나 branding 설정을 바꾼 뒤 기존 UI를 한 번에 재시작하려면 다음을 사용합니다.

```bash
./launch_ui.sh --restart
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

1. `Reference 폴더 탐색`에서 폴더 버튼을 눌러 이동하고, FASTA가 들어 있는 폴더에서
   `이 폴더 사용`을 누릅니다. 탐색 시작 위치와 최상위 범위는 `/data`입니다.
2. Reference는 파일명 기준 자연 정렬되며(예: `01, 02, ..., 10`),
   필요할 때만 접힌 `Reference 확인`에서 파일명과 길이를 확인합니다.
3. `ONT 결과 폴더 탐색`에서 sample 하위 폴더가 들어 있는 상위 폴더를 선택합니다.
   표준 ONT run 폴더의 `fastq_pass`만 자동 검색하며 `fastq_fail`, `other_reports`,
   `pod5_*`는 제외합니다. FASTQ/FASTQ.GZ 원본은 복사하지 않고 직접 연결합니다.
4. `Reference별 ONT sample 선택`에서 분석할 조합을 직접 지정합니다. 이름이나 순서를
   기준으로 자동 배정하지 않습니다. Reference마다 sample 수를 다르게 선택할 수 있고,
   같은 sample을 여러 Reference와 비교할 수도 있습니다. 선택하지 않은 Reference와
   ONT sample은 이번 분석에서 제외됩니다.
5. 왼쪽 하단의 `Analysis settings` 버튼을 누르면 중앙에 나타나는 설정 화면에서
   CPU, 병렬 sample 수, read filter 및 variant review 기준을 저장할 수 있습니다.
6. 중복·누락 경고가 없는지 확인한 뒤 `Batch analysis 실행`을 누릅니다.

왼쪽 sidebar의 작은 `분석 도구` checklist에서 `minimap2`, `samtools`, `bcftools`,
`gzip`, `bash`가 모두 ✅인지 확인할 수 있습니다. `NanoFilt`는 optional QC이므로
설치되지 않으면 ⚪로 표시됩니다. Sidebar를 접어도 화면 왼쪽 위의 화살표로 다시
열 수 있습니다.

각 입력 항목의 `?`에 커서를 올리면 설정 의미와 권장 사용법을 확인할 수 있습니다.

Offline 환경에서 `pyarrow`와 `arrow-cpp`가 호환되지 않으면 UI는 자동으로 Markdown
table을 사용합니다. Pyarrow가 정상 import되는 환경에서는 interactive table과 depth
chart가 자동 활성화되며, 어느 경우에도 variant calling 결과에는 차이가 없습니다.

결과는 reference별로 연결된 sample이 묶여 표시되며 `CLEAN`,
`VARIANT DETECTED`, `REVIEW`, `ERROR` 상태를 제공합니다. 전체 결과 요약은 CSV로
내려받을 수 있고 BAM, VCF, consensus FASTA 및 보고서는 서버 결과 폴더에 남습니다.
실행 command, 상세 log, 서버 결과 경로는 UI에 표시되지 않습니다. 관리자 진단을 위한
log 파일은 서버 실행 폴더에 계속 보존됩니다.

Sample ID는 `barcode13` 같은 번호, ONT sample alias 폴더명 또는 FASTQ 파일명에서
자동으로 가져옵니다. 따라서 barcode라는 이름을 반드시 사용할 필요가 없습니다.

> 기본 탐색 범위는 `/data`입니다. 다른 범위가 꼭 필요한 경우에만 UI 실행 전에
> `ONT_SERVER_ROOT` 환경변수를 지정합니다.

## 4. 결과 위치

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

## 5. 정확도 해석

- SNP는 이 UI에서 1-bp substitution/point mutation을 의미합니다.
- `PASS`는 현재 선택한 QUAL, DP, allele-fraction, edge 기준을 만족합니다.
- `REVIEW`는 오류 확정이 아니라 추가 확인이 필요한 call입니다.
- ONT indel은 homopolymer에서 false positive가 증가할 수 있으므로 BAM을 IGV로
  열어 read-level support를 함께 확인하는 것이 좋습니다.
- plasmid/vector 단일 clone의 진짜 변이는 보통 높은 allele fraction을 보입니다.

## 6. 성능 최적화

- Pipeline command는 브라우저에 실시간으로 반복 출력하지 않고 `pipeline.log`에 저장합니다.
- 공유 filesystem에서 불필요한 source scan이 발생하지 않도록 Streamlit file watcher를 끕니다.
- 진행 화면은 완료된 sample 수만 갱신합니다.
- Sample ID 편집표와 reference별 상세 결과는 사용자가 요청할 때만 생성합니다.
- Analysis log는 사용자가 선택할 때 최근 300줄만 표시합니다.

## 7. 테스트

```bash
python -m unittest discover -s tests -v
bash -n run_pipeline.sh run_ui.sh launch_ui.sh install_linux_launcher.sh
```
