# ONT Plasmid Analyzer example data

이 데이터는 UI와 pipeline 실행 확인을 위한 synthetic sequence입니다.

1. UI의 Reference 탐색기에서 `references` 폴더를 선택합니다.
2. ONT 결과 탐색기에서 `demo_ont_run` 폴더를 선택합니다.
3. `expected_variants.csv`를 보고 각 Reference에 ONT sample 3개를 직접 선택합니다.
4. `Batch analysis 실행`을 누릅니다.
5. 분석 결과를 `expected_variants.csv`와 비교합니다.

`demo_ont_run`은 실제 ONT output과 비슷하게 `fastq_pass`, `fastq_fail`,
`other_reports`, `pod5_fail`, `pod5_pass`, `pod5_skip` 폴더로 구성됩니다.
분석기는 `fastq_pass`의 FASTQ.GZ만 읽습니다. Sample 이름에는 예상 variant를
표시하지 않았으며 정답은 `expected_variants.csv`에서만 확인할 수 있습니다.

기본 Analysis settings에서 각 sample은 12 reads이므로 minimum depth 10 조건을 통과합니다.
실제 생물학적 데이터가 아닌 소프트웨어 테스트 전용 데이터입니다.
