# ONT Plasmid Analyzer example data

이 데이터는 UI와 pipeline 실행 확인을 위한 synthetic sequence입니다.

1. `references` 안의 FASTA 5개를 한 번에 업로드합니다.
2. 각 reference 영역에 같은 이름의 `samples/demo_plasmid_XX` 폴더를 드래그합니다.
3. Reference별로 Control, SNP, Insertion 또는 Deletion sample 3개가 감지되는지 확인합니다.
4. `Batch analysis 실행`을 누릅니다.
5. 결과를 `expected_variants.csv`와 비교합니다.

기본 Analysis settings에서 각 sample은 12 reads이므로 minimum depth 10 조건을 통과합니다.
실제 생물학적 데이터가 아닌 소프트웨어 테스트 전용 데이터입니다.
