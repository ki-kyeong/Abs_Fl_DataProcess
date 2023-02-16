# Abs_Fl_DataProcess
 
Functions for processing absorption and Fluorenscence time trace signal data of Li or MgF.

## Code Descriptions
>### Data processing codes
>>- #### data2spec.m
>>	초기 버전의 코드. csv file들을 불러와 주파수별로 처리해준다. 현재는 사용하지 않음.

>>- #### readLiSpecData.m
>>	Li experiment data용 코드. 7Li D2 F=2 transition frequency 417.80987 THz를 기준을 계산한다.
   
>>- #### readMgFData.m
>>	MgF experiment data용 코드. 여러 MgF transition line으 선택할 수 있게 만들 예정
   
>### Fitting codes
>>- #### GF_Li.m
>> Li absorption spectrum fitting code. single/double gaussian 둘다 가능.

>### Plot codes
>>- #### plotspec.m
>>	data structure를 받아 abs/fl spectrum을 그려줌
>>- #### plot2Dtimedet.m
>>	data structure를 받아 velocity vs time 2D image 생성.
>>- #### plotfigureset.m
>> data structure를 받아 abs/fl spectrum과 velocity vs time 2D image들을 한번에 생성
