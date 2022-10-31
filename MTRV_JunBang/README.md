# Multi Target Rendezvous Two-Phase FrameWork Research & Coded by Jun Bang

Jun Bang Original Research 중 source code만 저장

1. Case Study GTOC & IRIDIUM-Long의 경우 J2 Perturbation 을 활용한 Rendezvous (MATLAB)
2. Case Study IRIDIUM-Short 의 경우 Lambert Rendezvous를 활용함. (C)


Build 생성하는 방법
1. 새롭게 생성하고자 하는 이름의 프로젝트 생성 (콘솔 응용 프로그램 / 솔루션용 디렉터리 만들기 해제/ 다음 선택 후 빈 프로젝트 생성)
2. 기존 원본 폴더와 같은 형태로 구성
3. 프로젝트 소스파일에 해당 파일 추가 (원본 Source 폴더 내 확장자 .c)
4. 프로젝트 헤더파일에 해당 파일 추가 (원본 Source 폴더 내 확장자 .h)
5. SNOPT 라이브러리를 생성한 "프로젝트" 하위에 추가
6. 솔루션 우클릭 -> 구성관리자 -> Release, 새로만들기 x64
7. 프로젝트 우클릭 -> 속성 -> 구성속성 -> 링커 -> 입력 -> 특정 기본라이브러리 무시에
"libcmt.lib" 입력


C 코드의 VS2013 버전으로 Compile 해야함. 이후 버전 빌드 깨짐. 최신화 -> MTRV_ElementarySolGenaration
