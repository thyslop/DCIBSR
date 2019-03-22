*****************************************************************************************************************
* STUDY:	The impact of receipt and timing of adjuvant systemic therapy on outcomes following diagnosis of
*         node-positive invasive lobular carcinoma
* FILE NAME: inv_lobular_data.sas
* PURPOSE:	To determine if patients with node-positive invasive lobular disease who receive neoadjuvant
            chemotherapy have higher likelihood of death than patients with node-positive invasive lobular
            disease who do not receive neoadjuvant chemotherapy.
* UPDATE:	N/A
* PIs:	Lola Fayanju, MD
* DEPARTMENT:	Breast
* DATE CREATED:	29JAN2018
* DATE UPDATED: 11APR2018: updated to clean up some code from abstract and output new working dataset
*               03MAY2018: updated to remove pts with metastatic disease
*               14MAY2018: updated to exclude cT0, TX pts in primary analysis, updated NACT grouping
*               24MAY2018: updated exclusion criteria
*               08JUN2018: updated to remove patients missing cT from analysis set
*               05OCT2018: updated facility type and location derivations
* CREATED BY:	Hannah Vernia (HV)
* DATA:
*****************************************************************************************************************;
/****************************************************/
/* Set up for data/formats/directory                */
/****************************************************/
%let DIR=H:\Breast\Inv_Lobular\;

libname ncdb "H:\Breast\NCDB\";
libname proj "&dir.";
libname mac "H:\Macros";
options fmtsearch=(proj ncdb mac) minoperator;

/****************************************************/
/* Import data                                      */
/****************************************************/
data bcdata;
	set ncdb.bc2014;
run;
/*N=2246363*/

/****************************************************/
/* BEGIN I/E and population definitions             */
/****************************************************/
*Subsetting based on inclusion criteria - patients with node-positive lobular carcinoma;
data bcdata_ILC;
	set bcdata;
  if HISTOLOGY in (8520 8522 8524)
     and BEHAVIOR = 3
     and TNM_CLIN_N not in ("" "0" "X");
run;
/*N=34568*/

/****************************************************/
/* BEGIN populations using only ILC patients        */
/****************************************************/
*checking dropouts by each exclusion criterion;
proc sql;
  select count(distinct puf_case_id) as mets from bcdata_ILC where TNM_CLIN_M ne "0";
  select count(distinct puf_case_id) as tum from bcdata_ILC where TNM_CLIN_T in ("0" "X" "" "IS");
  select count(distinct puf_case_id) as tum4d from bcdata_ILC where index(TNM_CLIN_T,"4D")>0;
  select count(distinct puf_case_id) as chemo from bcdata_ILC where RX_SUMM_CHEMO not in (1 2 3);
  select count(distinct puf_case_id) as surg from bcdata_ILC where RX_SUMM_SURG_PRIM_SITE in (. 00 19 90 99);
  select count(distinct puf_case_id) as surv from bcdata_ILC where missing(puf_vital_status) or missing(dx_lastcontact_death_months);
  select count(distinct puf_case_id) as trtseq from bcdata_ILC where RX_SUMM_CHEMO in (1 2 3) and nmiss(DX_DEFSURG_STARTED_DAYS,DX_CHEMO_STARTED_DAYS)>0 and RX_SUMM_SYSTEMIC_SUR_SEQ in(. 0 9);
  select count(distinct puf_case_id) as naet from bcdata_ILC where .<DX_HORMONE_STARTED_DAYS<=DX_DEFSURG_STARTED_DAYS;
quit;

*Exclusion criteria - primary endpoint;
*5/3/18: updated to exclude mets for all populations;
*5/24/18: updated to exclude cT4d, neoadjuvant endocrine therapy;
data bcdata_ILC2;
  set bcdata_ILC;
  if TNM_CLIN_M ne "0" then delete; *remove pts with metastatic cancer;
  if TNM_CLIN_T in ("0" "X" "" "IS") then delete; *remove pts with cT0, TX, Tis, missing cT;
  if index(TNM_CLIN_T,"4D")>0 then delete; *remove pts with cT4d;
  if RX_SUMM_CHEMO in (0 1 2 3 82 85 86 87); *Keep only pts with known chemo status;
  if RX_SUMM_SURG_PRIM_SITE in (. 00 19 90 99) then delete; *Remove pts who did not undergo lumpectomy or mastectomy;
  if missing(PUF_VITAL_STATUS) or missing(DX_LASTCONTACT_DEATH_MONTHS) then delete; *Remove pts with missing survival;
  if RX_SUMM_CHEMO in (1 2 3) and nmiss(DX_DEFSURG_STARTED_DAYS, DX_CHEMO_STARTED_DAYS)>0 and RX_SUMM_SYSTEMIC_SUR_SEQ in(. 0 9) then delete; *Remove pts with missing info on treatment sequence;
  if .<DX_HORMONE_STARTED_DAYS<=DX_DEFSURG_STARTED_DAYS then delete; *remove pts with neoadjuvant endocrine therapy;
run;
proc sort presorted; by PUF_CASE_ID; run;
/*N=19924*/

*Analysis groups - primary endpoint: cN1+ ILC: NACT vs. no NACT;
data bcdata_prim;
  set bcdata_ILC2 (where=(RX_SUMM_CHEMO in (1 2 3)));
  length TRT01A $200;
  if .<DX_CHEMO_STARTED_DAYS<=DX_DEFSURG_STARTED_DAYS or (RX_SUMM_SYSTEMIC_SUR_SEQ in (2 4) and (nmiss(DX_CHEMO_STARTED_DAYS,DX_DEFSURG_STARTED_DAYS)>0 or DX_CHEMO_STARTED_DAYS<=DX_DEFSURG_STARTED_DAYS)) then do;
    TRT01A = "NACT";
    TRT01AN = 1;
  end;
  else do;
    TRT01A = "ACT";
    TRT01AN = 2;
  end;
run;
proc freq data=bcdata_prim (where=(not missing(trt01a))); tables trt01a; run;
/*N=15573*/

*Exclusion criteria and analysis groups - secondary endpoint #1: cN1-3 ILC: <50 vs >=50;
*NOTE: this may need to be updated to exclude patients who didn't receive chemo;
data bcdata_sec1;
  set bcdata_ILC2 (where=(not missing(AGE)));
  length TRT02A $200;
  if . < AGE < 50 then do;
    TRT02A = "Age < 50";
    TRT02AN = 1;
  end;
  else if AGE >= 50 then do;
    TRT02A = "Age >= 50";
    TRT02AN = 2;
  end;
run;
proc freq data=bcdata_sec1 (where=(not missing(trt02a))); tables trt02a; run;
/*N=19924*/

*Exclusion criteria and analysis groups - secondary endpoint #2: cN1-3 ILC: chemo vs. endo vs. chemo+endo;
data bcdata_sec2;
  set bcdata_ILC2 (where=(RX_SUMM_HORMONE=1 or RX_SUMM_CHEMO in (1 2 3)));
  length TRT03A $200;
  if RX_SUMM_HORMONE ne 1 and RX_SUMM_CHEMO in (1 2 3) then do;
    TRT03A = "Chemotherapy only";
    TRT03AN = 1;
  end;
  else if RX_SUMM_HORMONE = 1 and RX_SUMM_CHEMO not in (1 2 3) then do;
    TRT03A = "Endocrine therapy only";
    TRT03AN = 2;
  end;
  else if RX_SUMM_HORMONE = 1 and RX_SUMM_CHEMO in (1 2 3) then do;
    TRT03A = "Chemotherapy and Endocrine therapy";
    TRT03AN = 3;
  end;
run;
proc freq data=bcdata_sec2 (where=(not missing(trt03a))); tables trt03a; run;
/*N=18598*/

/****************************************************/
/* END populations using only ILC patients          */
/* BEGIN populations using ILC and IDC patients     */
/****************************************************/
*Inclusion criteria - IDC patients;
*5/3/18: updated to exclude mets for all populations;
data bcdata_IDC;
	set bcdata;
  if HISTOLOGY in (8500 8501 8502 8503 8504 8507 8508 8509 8022 8035 8201 8211 8480 8540)
     and (CS_SITESPECIFIC_FACTOR_1=10 or CS_SITESPECIFIC_FACTOR_2=10 or CS_SITESPECIFIC_FACTOR_16 in (10 11 100 101 110 111))
     and BEHAVIOR = 3
     and TNM_CLIN_N not in ("" "0" "X");
run;
/*N=135499*/

*Exclusion criteria - IDC patients;
data bcdata_IDC2;
  set bcdata_IDC;
  if TNM_CLIN_M = "0" then delete; *drop pts with mets;
  if TNM_CLIN_T in ("0" "X" "" "IS") then delete; *drop pts with cT0, TX;
  if index(TNM_CLIN_T,"4D")>0 then delete; *remove pts with cT4d;
  if RX_SUMM_CHEMO in (0 1 2 3 82 85 86 87); *Keep only pts with known chemo status;
  if RX_SUMM_SURG_PRIM_SITE in (. 00 19 90 99) then delete; *Remove pts who did not undergo lumpectomy or mastectomy;
  if missing(PUF_VITAL_STATUS) or missing(DX_LASTCONTACT_DEATH_MONTHS) then delete; *Remove pts with missing survival;
  *if RX_SUMM_CHEMO in (1 2 3) and (missing(DX_DEFSURG_STARTED_DAYS) or missing(DX_CHEMO_STARTED_DAYS)) then delete; *Remove pts with missing info on treatment sequence;
  if .<DX_HORMONE_STARTED_DAYS<=DX_DEFSURG_STARTED_DAYS then delete; *remove pts with neoadjuvant endocrine therapy;
run;
proc sort presorted; by PUF_CASE_ID; run;
/*N=6201*/

*Analysis groups - secondary endpoint #3: Compare cN -> pN for cN1-3 ILC to cN -> pN for cN1-3 IDC;
*NOTE: if analysis will be done on these populations, please confirm I/E criteria;
data bcdata_sec3;
  set bcdata_ILC2 (in=inILC) bcdata_IDC2 (in=inIDC);

  length TRT04A $200;
  if inILC then do;
    TRT04A = "Invasive lobular carcinoma";
    TRT04AN = 1;
  end;
  else if inIDC then do;
    TRT04A = "Invasive ductal carcinoma";
    TRT04AN = 2;
  end;

  if missing(TNM_CLIN_N) or missing(TNM_PATH_N) then delete; *Remove pts where stage change cannot be assessed;
run;
proc sort presorted; by PUF_CASE_ID; run;
proc freq data=bcdata_sec3 (where=(not missing(trt04a))); tables trt04a; run;
/*N=25648*/

*Analysis groups - secondary endpoint #4: pts <50: cN1-3 ILC vs cN1-3 IDC;
data bcdata_sec4;
  set bcdata_ILC2 (in=inILC) bcdata_IDC2 (in=inIDC);

  length TRT05A $200;
  if inILC then do;
    TRT05A = "Invasive lobular carcinoma";
    TRT05AN = 1;
  end;
  else if inIDC then do;
    TRT05A = "Invasive ductal carcinoma";
    TRT05AN = 2;
  end;

  if AGE >=50 then delete; *Only keep pts < 50yr for this analysis;
run;
proc sort presorted; by PUF_CASE_ID; run;
proc freq data=bcdata_sec4 (where=(not missing(trt05a))); tables trt05a; run;
/*N=7006*/

/****************************************************/
/* END populations using only ILC patients          */
/* BEGIN creation of permanent dataset              */
/****************************************************/
data bcdata1;
  merge bcdata_ILC2
        bcdata_IDC2
        bcdata_prim (keep=PUF_CASE_ID TRT01A TRT01AN)
        bcdata_sec1 (keep=PUF_CASE_ID TRT02A TRT02AN)
        bcdata_sec2 (keep=PUF_CASE_ID TRT03A TRT03AN)
        bcdata_sec3 (keep=PUF_CASE_ID TRT04A TRT04AN)
        bcdata_sec4 (keep=PUF_CASE_ID TRT05A TRT05AN);
  by PUF_CASE_ID;
  if nmiss(TRT01AN, TRT02AN, TRT03AN, TRT04AN, TRT05AN)<5;
run;
/*N=25953*/

/****************************************************/
/* END I/E and population definitions               */
/* BEGIN analysis programming - variable definitions*/
/****************************************************/
*Per NCDB guidelines, hospital location and type are missing for patients < 40, so creating a dataset that has
	hosptial ID, type, location for merge onto dataset;
proc sql;
	create table hosps as select count(unique(puf_case_id)) as numpts, puf_facility_id, FACILITY_TYPE_CD, FACILITY_LOCATION_CD 
		from bcdata group by	puf_facility_id, FACILITY_TYPE_CD, FACILITY_LOCATION_CD;
quit;
data hosps;
	set hosps;
	if FACILITY_TYPE_CD=FACILITY_LOCATION_CD=. then delete;
run;
proc sort data=hosps;
	by puf_facility_id descending numpts;
run;

data hosps2;
  set hosps;
  by puf_facility_id descending numpts;
  if first.puf_facility_id;
run;

*merging these facility types and locations onto the data;
proc sort data=bcdata1; by puf_facility_id; run;
data bcdata2;
	merge bcdata1 hosps2(rename=(FACILITY_TYPE_CD=FACILITY_TYPE_CD2 FACILITY_LOCATION_CD=FACILITY_LOCATION_CD2));
	by puf_facility_id;
run;

*Redefining variables of interest;
data bcdata3;
	length sex2 loc subtype $10. age2 inc cd_score tum fac ins ethnicity laterality2 comb_erprher2 $25. race2 edu radtype $50.;
	set bcdata2;

	if CS_SITESPECIFIC_FACTOR_1=10 or CS_SITESPECIFIC_FACTOR_16 in (100 101 110 111) then ER="ER+";
		else if CS_SITESPECIFIC_FACTOR_1=20 or CS_SITESPECIFIC_FACTOR_16 in (0 1 10 11) then ER="ER-";
	if CS_SITESPECIFIC_FACTOR_2=10 or CS_SITESPECIFIC_FACTOR_16 in (10 11 110 111) then PR="PR+";
		else if CS_SITESPECIFIC_FACTOR_2=20 or CS_SITESPECIFIC_FACTOR_16 in (0 1 100 101) then PR="PR-";
	if CS_SITESPECIFIC_FACTOR_15=10 or CS_SITESPECIFIC_FACTOR_16 in (1 11 101 111) then HER2="HER2+";
		else if CS_SITESPECIFIC_FACTOR_15=20 or CS_SITESPECIFIC_FACTOR_16 in (0 10 100 110) then HER2="HER2-";
	if er="" or pr="" or her2="" then comb_erprher2="";
		else comb_erprher2=strip(ER)||"/"||strip(PR)||"/"||strip(HER2);

  if comb_erprher2="ER-/PR-/HER2-" then subtype="TNBC";
		else if comb_erprher2 in ("ER+/PR+/HER2-" "ER+/PR-/HER2-" "ER-/PR+/HER2-") then subtype="Luminal";
		else if comb_erprher2 in ("ER+/PR+/HER2+" "ER+/PR-/HER2+" "ER-/PR+/HER2+" "ER-/PR-/HER2+") then subtype="HER2+";

	if SEX=1 then sex2="Male";
		else if SEX=2 then sex2="Female";

	if race=1 then race2="White";
		else if race=2 then race2="Black";
		else if race=99 then race2="";
		else if race in (3 4 5 6 8 10 11 12 13 14 15 16 17 96 7 20 21 22 25 26 27 28 30 31 32 97 98) then race2="Other";

	if SPANISH_HISPANIC_ORIGIN=0 then ethnicity="Non-Hispanic";
		else if SPANISH_HISPANIC_ORIGIN in (1 2 3 4 5 6 7 8) then ethnicity="Hispanic";

	if MED_INC_QUAR_00 in (1 2) then inc="<$35,000";
		else if MED_INC_QUAR_00 in (3 4) then inc="^{unicode '2265'x}$35,000";

	if NO_HSD_QUAR_00 in (1 2) then edu="^{unicode '2264'x}80% High School Graduation Rate";
		else if NO_HSD_QUAR_00 in (3 4) then edu=">80% High School Graduation Rate";

	if INSURANCE_STATUS=0 then ins="Not Insured";
		else if INSURANCE_STATUS in (2 3 4) then ins="Government";
		else if INSURANCE_STATUS=1 then ins="Private";

	if CDCC_TOTAL=0 then cd_score="0";
		else if CDCC_TOTAL=1 then cd_score="1";
		else if CDCC_TOTAL>=2 then cd_score="^{unicode '2265'x}2";

	if FACILITY_TYPE_CD2=1 then fac="Community";
		else if FACILITY_TYPE_CD2=2 then fac="Comprehensive";
		else if FACILITY_TYPE_CD2=3 then fac="Academic";
		else if FACILITY_TYPE_CD2=4 then fac="Integrated Network";
		else if FACILITY_TYPE_CD2=9 then fac="Other";

	if FACILITY_LOCATION_CD2 in (1 2) then loc="Northeast";
		else if FACILITY_LOCATION_CD2 in (3 5 7) then loc="South";
		else if FACILITY_LOCATION_CD2 in (4 6) then loc="Midwest";
		else if FACILITY_LOCATION_CD2 in (8 9) then loc="West";

	if 0<tumor_size<=10 or tumor_size=991 then tum="<1 cm";
		else if 11<=tumor_size<=20 or tumor_size=992 then tum=">1 to 2 cm";
		else if 21<=tumor_size<=40 or tumor_size in (993 994) then tum=">2 to 4 cm";
		else if 41<=tumor_size<=990 or tumor_size=995 then tum=">4 cm";

	if 0<tumor_size<991 then tum2=tumor_size;
	if not missing(tum2) then tum_cm=tum2/10;

	distance=crowfly;

	if LATERALITY in (1 2 3) then laterality2="Unilateral";
		else if LATERALITY in (4 5) then laterality2="Bilateral";

	if RX_SUMM_RADIATION in (1 2 3 4 5) then rad="Yes";
		else if RX_SUMM_RADIATION=0 then rad="No";

	if RX_SUMM_CHEMO in (1 2 3) then chemo="Yes";
		else if RX_SUMM_CHEMO in (0 82 85 86 87) then chemo="No";

	os_mos=DX_LASTCONTACT_DEATH_MONTHS;

	if puf_vital_status=1 then dead=0;
		else if puf_vital_status=0 then dead=1;

  if index(TNM_CLIN_T,"0")>0 then clin_t="0";
    else if index(TNM_CLIN_T,"1")>0 then clin_t="1";
		else if index(TNM_CLIN_T,"2")>0 then clin_t="2";
    else if index(TNM_CLIN_T,"3")>0 then clin_t="3";
    else if index(TNM_CLIN_T,"4")>0 then clin_t="4";
    else if index(TNM_CLIN_T,"X")>0 then clin_t="X";

  if index(TNM_CLIN_N,"0")>0 then clin_n="0";
    else if index(TNM_CLIN_N,"1")>0 then clin_n="1";
		else if index(TNM_CLIN_N,"2")>0 then clin_n="2";
    else if index(TNM_CLIN_N,"3")>0 then clin_n="3";
    else if index(TNM_CLIN_N,"X")>0 then clin_n="X";

  if index(TNM_CLIN_M,"0")>0 then clin_m="0";
    else if index(TNM_CLIN_M,"1")>0 then clin_m="1";
		else if index(TNM_CLIN_M,"2")>0 then clin_m="2";
    else if index(TNM_CLIN_M,"3")>0 then clin_m="3";
    else if index(TNM_CLIN_M,"X")>0 then clin_m="X";

  *Adding pathological T and N;
    if index(TNM_PATH_T,"0")>0 then path_t="0";
    else if index(TNM_PATH_T,"1")>0 then path_t="1";
		else if index(TNM_PATH_T,"2")>0 then path_t="2";
    else if index(TNM_PATH_T,"3")>0 then path_t="3";
    else if index(TNM_PATH_T,"4")>0 then path_t="4";
    else if index(TNM_PATH_T,"X")>0 then path_t="X";

  if index(TNM_PATH_N,"0")>0 then path_n="0";
    else if index(TNM_PATH_N,"1")>0 then path_n="1";
		else if index(TNM_PATH_N,"2")>0 then path_n="2";
    else if index(TNM_PATH_N,"3")>0 then path_n="3";
    else if index(TNM_PATH_N,"X")>0 then path_n="X";

  *Adding stages II/III for sensitivity analysis (based on AJCC 7th ed);
  if index(TNM_CLIN_T,"1")>0 then do;
    if index(TNM_CLIN_N,"1")>0 then stage=2;
    else if index(TNM_CLIN_N,"2")>0 or index(TNM_CLIN_N,"3")>0 then stage=3;
  end;
  else if index(TNM_CLIN_T,"2")>0 then do;
    if index(TNM_CLIN_N,"1")>0 then stage=2;
    else if index(TNM_CLIN_N,"2")>0 or index(TNM_CLIN_N,"3")>0 then stage=3;
  end;
  else if index(TNM_CLIN_T,"3")>0 then do;
    if index(TNM_CLIN_N,"1")>0 then stage=2;
    else if index(TNM_CLIN_N,"2")>0 or index(TNM_CLIN_N,"3")>0 then stage=3;
  end;
  else if index(TNM_CLIN_T,"4")>0 then stage=3;

	if .<age<50 then age2="<50";
	else if age>=50 then age2="^{unicode '2265'x}50";

	if grade=1 then grade2="1";
		else if grade=2 then grade2="2";
		else if grade in (3 4) then grade2="3";	*include 3 and 4 because breast cancer can't have grade 4;
		else if grade=9 then grade2="";

	if 0<=CS_SITESPECIFIC_FACTOR_3<=90 then num_pos_axil_nodes=CS_SITESPECIFIC_FACTOR_3;
		else if CS_SITESPECIFIC_FACTOR_3 in (. 95 97 98 99) and 0<=REGIONAL_NODES_POSITIVE<=90 then num_pos_axil_nodes=REGIONAL_NODES_POSITIVE;

	if 1<=REGIONAL_NODES_EXAMINED<=90 then num_nodes_exam=REGIONAL_NODES_EXAMINED;

	if RX_SUMM_HORMONE=1 then hormone=1;
		else if RX_SUMM_HORMONE in (0 82 85 86 87) then hormone=0;

  if not missing(year_of_diagnosis) then year_of_diagnosis_c=put(year_of_diagnosis,4.);

  if RX_SUMM_SURG_PRIM_SITE in (20 21 22 23 24) then surgtype = "Lumpectomy";
  else surgtype = "Mastectomy";

  if surgtype = "Lumpectomy" then do;
    if DX_RAD_STARTED_DAYS > DX_DEFSURG_STARTED_DAYS then radtype = "Post-Lumpectomy RT";
    else if rad="Yes" then radtype = "Other RT";
    else radtype = "No RT";    
  end;
  else do;
    if DX_RAD_STARTED_DAYS > DX_DEFSURG_STARTED_DAYS then radtype = "Post-Mastectomy RT";
    else if rad="Yes" then radtype = "Other RT";
    else radtype = "No RT";    
  end;
run;

*Checking that these variables were made correctly;
proc freq data=bcdata3;
	table tnm_clin_t*clin_t tnm_clin_m*clin_m tnm_clin_n*clin_n grade*grade2 RX_SUMM_HORMONE*hormone puf_vital_status*dead rx_summ_radiation*rad rx_summ_chemo*chemo 
		laterality*laterality2 comb_erprher2*subtype age*age2 rx_summ_radiation radtype / missing;
run;

/****************************************************/
/* Create permanent output dataset                  */
/****************************************************/
data proj.analysisdata&sysdate9.;
	set bcdata3;
run;
