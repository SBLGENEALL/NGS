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

1. Windows Explorer의 Reference FASTA는 drag-and-drop 영역에 올립니다. MobaXterm
   SFTP에서 브라우저로 직접 드래그되지 않는 환경에서는 아래의 접힌
   `Linux server에 있는 Reference 불러오기`를 열고 `pwd`로 확인한 파일 또는 폴더의
   절대경로를 입력합니다. PC 파일과 server 파일을 함께 사용할 수도 있습니다.
2. Reference는 파일명 기준 자연 정렬되며(예: `01, 02, ..., 10`),
   `Reference file check`에서 정렬된 파일명과 길이를 확인합니다.
3. Reference별로 Windows 업로드 또는 Linux server 경로를 선택합니다. Server 경로는
   FASTQ 한 개 또는 sample 하위 폴더를 포함한 상위 폴더를 지정할 수 있습니다.
   Server FASTQ는 복사하지 않고 원본을 직접 연결하므로 대용량 upload 시간이 없습니다.
4. 필요하면 `Sample ID 직접 수정`을 켜서 자동 인식된 ONT sample name을 수정합니다.
   여러 FASTQ chunk에 같은 Sample ID를 입력하면 하나의 sample로 합쳐 분석합니다.
5. `분석 전 최종 확인`에서 reference별 Sample ID와 FASTQ 개수를 확인합니다.
6. 왼쪽의 `Analysis settings` 버튼을 누르면 중앙에 나타나는 설정 화면에서
   CPU, 병렬 sample 수, read filter 및 variant review 기준을 저장할 수 있습니다.
7. 중복·누락 경고가 없는지 확인한 뒤 `Batch analysis 실행`을 누릅니다.

각 입력 항목의 `?`에 커서를 올리면 설정 의미와 권장 사용법을 확인할 수 있습니다.

결과는 reference별로 연결된 sample이 묶여 표시되며 `CLEAN`,
`VARIANT DETECTED`, `REVIEW`, `ERROR` 상태를 제공합니다. 전체 결과 요약은 CSV로
내려받을 수 있고 BAM, VCF, consensus FASTA 및 보고서는 서버 결과 폴더에 남습니다.
실행 command와 상세 출력은 기본적으로 접힌 `Analysis log`에서 필요할 때만 확인합니다.

Sample ID는 `barcode13` 같은 번호, ONT sample alias 폴더명 또는 FASTQ 파일명에서
자동으로 가져옵니다. 따라서 barcode라는 이름을 반드시 사용할 필요가 없습니다.

> 디렉터리 업로드를 위해 Streamlit 1.57 이상이 필요합니다. FASTQ는 브라우저를
> 통해 서버로 전송되므로 분석 중 브라우저 탭을 닫지 마세요.

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
