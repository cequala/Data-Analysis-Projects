/* ----------------------------------------
Code exported from SAS Enterprise Guide
DATE: 2025³â 8¿ù 12ÀÏ     TIME: 20:47:15
PROJECT: GangMira
PROJECT PATH: C:\Users\Admin\Desktop\sas_code\GangMira.egp
---------------------------------------- */

/* ---------------------------------- */
/* MACRO: enterpriseguide             */
/* PURPOSE: define a macro variable   */
/*   that contains the file system    */
/*   path of the WORK library on the  */
/*   server.  Note that different     */
/*   logic is needed depending on the */
/*   server type.                     */
/* ---------------------------------- */
%macro enterpriseguide;
%global sasworklocation;
%local tempdsn unique_dsn path;

%if &sysscp=OS %then %do; /* MVS Server */
	%if %sysfunc(getoption(filesystem))=MVS %then %do;
        /* By default, physical file name will be considered a classic MVS data set. */
	    /* Construct dsn that will be unique for each concurrent session under a particular account: */
		filename egtemp '&egtemp' disp=(new,delete); /* create a temporary data set */
 		%let tempdsn=%sysfunc(pathname(egtemp)); /* get dsn */
		filename egtemp clear; /* get rid of data set - we only wanted its name */
		%let unique_dsn=".EGTEMP.%substr(&tempdsn, 1, 16).PDSE"; 
		filename egtmpdir &unique_dsn
			disp=(new,delete,delete) space=(cyl,(5,5,50))
			dsorg=po dsntype=library recfm=vb
			lrecl=8000 blksize=8004 ;
		options fileext=ignore ;
	%end; 
 	%else %do; 
        /* 
		By default, physical file name will be considered an HFS 
		(hierarchical file system) file. 
		*/
		%if "%sysfunc(getoption(filetempdir))"="" %then %do;
			filename egtmpdir '/tmp';
		%end;
		%else %do;
			filename egtmpdir "%sysfunc(getoption(filetempdir))";
		%end;
	%end; 
	%let path=%sysfunc(pathname(egtmpdir));
    %let sasworklocation=%sysfunc(quote(&path));  
%end; /* MVS Server */
%else %do;
	%let sasworklocation = "%sysfunc(getoption(work))/";
%end;
%if &sysscp=VMS_AXP %then %do; /* Alpha VMS server */
	%let sasworklocation = "%sysfunc(getoption(work))";                         
%end;
%if &sysscp=CMS %then %do; 
	%let path = %sysfunc(getoption(work));                         
	%let sasworklocation = "%substr(&path, %index(&path,%str( )))";
%end;
%mend enterpriseguide;

%enterpriseguide


/* Conditionally delete set of tables or views, if they exists          */
/* If the member does not exist, then no action is performed   */
%macro _eg_conditional_dropds /parmbuff;
	
   	%local num;
   	%local stepneeded;
   	%local stepstarted;
   	%local dsname;
	%local name;

   	%let num=1;
	/* flags to determine whether a PROC SQL step is needed */
	/* or even started yet                                  */
	%let stepneeded=0;
	%let stepstarted=0;
   	%let dsname= %qscan(&syspbuff,&num,',()');
	%do %while(&dsname ne);	
		%let name = %sysfunc(left(&dsname));
		%if %qsysfunc(exist(&name)) %then %do;
			%let stepneeded=1;
			%if (&stepstarted eq 0) %then %do;
				proc sql;
				%let stepstarted=1;

			%end;
				drop table &name;
		%end;

		%if %sysfunc(exist(&name,view)) %then %do;
			%let stepneeded=1;
			%if (&stepstarted eq 0) %then %do;
				proc sql;
				%let stepstarted=1;
			%end;
				drop view &name;
		%end;
		%let num=%eval(&num+1);
      	%let dsname=%qscan(&syspbuff,&num,',()');
	%end;
	%if &stepstarted %then %do;
		quit;
	%end;
%mend _eg_conditional_dropds;


/* save the current settings of XPIXELS and YPIXELS */
/* so that they can be restored later               */
%macro _sas_pushchartsize(new_xsize, new_ysize);
	%global _savedxpixels _savedypixels;
	options nonotes;
	proc sql noprint;
	select setting into :_savedxpixels
	from sashelp.vgopt
	where optname eq "XPIXELS";
	select setting into :_savedypixels
	from sashelp.vgopt
	where optname eq "YPIXELS";
	quit;
	options notes;
	GOPTIONS XPIXELS=&new_xsize YPIXELS=&new_ysize;
%mend _sas_pushchartsize;

/* restore the previous values for XPIXELS and YPIXELS */
%macro _sas_popchartsize;
	%if %symexist(_savedxpixels) %then %do;
		GOPTIONS XPIXELS=&_savedxpixels YPIXELS=&_savedypixels;
		%symdel _savedxpixels / nowarn;
		%symdel _savedypixels / nowarn;
	%end;
%mend _sas_popchartsize;


%*--------------------------------------------------------------*
 * Tests the current version against a required version. A      *
 * negative result means that the SAS server version is less    *
 * than the version required.  A positive result means that     *
 * the SAS server version is greater than the version required. *
 * A result of zero indicates that the SAS server is exactly    *
 * the version required.                                        *
 *                                                              *
 * NOTE: The parameter maint is optional.                       *
 *--------------------------------------------------------------*;
%macro _SAS_VERCOMP(major, minor, maint);
    %_SAS_VERCOMP_FV(&major, &minor, &maint, &major, &minor, &maint)
%mend _SAS_VERCOMP;

%*--------------------------------------------------------------*
 * Tests the current version against either the required        *
 * foundation or Viya required version depending on whether the *
 * SYSVLONG version is a foundation or Viya one. A negative     *
 * result means that the SAS server version is less than the    *
 * version required.  A positive result means that the SAS      *
 * server version is greater than the version required. A       *
 * result of zero indicates that the SAS server is exactly the  *
 * version required.                                            *
 *                                                              *
 * NOTE: The *maint parameters are optional.                    *
 *--------------------------------------------------------------*;
%macro _SAS_VERCOMP_FV(fmajor, fminor, fmaint, vmajor, vminor, vmaint);
    %local major;
    %local minor;
    %local maint;
    %local CurMaj;
    %local CurMin;
    %local CurMnt;

    %* Pull the current version string apart.;
    %let CurMaj = %scan(&sysvlong, 1, %str(.));

    %* The Viya version number has a V on the front which means
       we need to adjust the Maint SCAN funtion index and also
       get the appropriate parameters for the major, minor, and
       maint values we need to check against (foundation or Viya);
    %if %eval(&CurMaj EQ V) %then
        %do;
		   %*   MM mm t           MM = Major version , mm = Minor version , t = Maint version ;
		   %* V.03.04M2P07112018 ;

            %let major = &vmajor;
            %let minor = &vminor;
            %let maint = &vmaint;
			%let CurMaj = %scan(&sysvlong, 2, %str(.));
			%* Index is purposely 2 because V is now one of the scan delimiters ;
			%let CurMin = %scan(&sysvlong, 2, %str(.ABCDEFGHIKLMNOPQRSTUVWXYZ));
			%let CurMnt = %scan(&sysvlong, 3, %str(.ABCDEFGHIKLMNOPQRSTUVWXYZ));
        %end;
    %else
        %do;
		    %* M mm    t           M = Major version , mm = Minor version , t = Maint version ;  
		    %* 9.01.02M0P11212005 ;

            %let major = &fmajor;
            %let minor = &fminor;
            %let maint = &fmaint;
			%let CurMin = %scan(&sysvlong, 2, %str(.));
			%let CurMnt = %scan(&sysvlong, 4, %str(.ABCDEFGHIKLMNOPQRSTUVWXYZ));
        %end;

    %* Now perform the version comparison.;
    %if %eval(&major NE &CurMaj) %then
        %eval(&CurMaj - &major);
    %else
        %if %eval(&minor NE &CurMin) %then
            %eval(&CurMin - &minor);
        %else
            %if "&maint" = "" %then
                %str(0);
            %else
                %eval(&CurMnt - &maint);
%mend _SAS_VERCOMP_FV;

%*--------------------------------------------------------------*
 * This macro calls _SAS_VERCONDCODE_FV() with the passed       *
 * version. If the current server version matches or is newer,  *
 * then the true code (tcode) is executed, else the false code  *
 * (fcode) is executed.                                         *
 * Example:                                                     *
 *  %let isV92 =                                                *
 *     %_SAS_VERCONDCODE(9,2,0,                                 *
 *         tcode=%nrstr(Yes),                                   *
 *         fcode=%nrstr(No))                                    *
 *--------------------------------------------------------------*;
%macro _SAS_VERCONDCODE( major, minor, maint, tcode=, fcode= );
    %_SAS_VERCONDCODE_FV( &major, &minor, &maint, &major, &minor, &maint, &tcode, fcode )
%mend _SAS_VERCONDCODE;

%*--------------------------------------------------------------*
 * This macro calls _SAS_VERCOMP_FV() with the passed versions. *
 * If the current server version matches or is newer, then the  *
 * true code (tcode) is executed, else the false code (fcode)   *
 * is executed.                                                 *
 * Example:                                                     *
 *  %let isV92 =                                                *
 *     %_SAS_VERCONDCODE_FV(9,2,0, 3,5,0                        *
 *         tcode=%nrstr(Yes),                                   *
 *         fcode=%nrstr(No))                                    *
 *--------------------------------------------------------------*;
%macro _SAS_VERCONDCODE_FV( fmajor, fminor, fmaint, vmajor, vminor, vmaint, tcode=, fcode= );
    %if %_SAS_VERCOMP_FV(&fmajor, &fminor, &fmaint, &vmajor, &vminor, &vmaint) >= 0 %then
        %do;
        &tcode
        %end;
    %else
        %do;
        &fcode
        %end;
%mend _SAS_VERCONDCODE_FV;

%*--------------------------------------------------------------*
 * Tests the current version to see if it is a Viya version     *
 * number.                                                      *
 * A result of 1 indicates that the SAS server is a Viya        *
 * server.                                                      *
 * A zero result indicates that the server version is not       *
 * that of a Viya server.                                       *
 *--------------------------------------------------------------*;
%macro _SAS_ISVIYA;
    %local Major;

    %* Get the major component of the current version string.;
    %let Major = %scan(&sysvlong, 1, %str(.));

    %* Check if it it V for Viya.;
    %if %eval(&Major EQ V) %then
        %str(1);
    %else
        %str(0);
%mend _SAS_ISVIYA;


ODS PROCTITLE;
OPTIONS DEV=SVG;
GOPTIONS XPIXELS=0 YPIXELS=0;
%macro HTML5AccessibleGraphSupported;
    %if %_SAS_VERCOMP_FV(9,4,4, 0,0,0) >= 0 %then ACCESSIBLE_GRAPH;
%mend;
FILENAME EGHTMLX TEMP;
ODS HTML5(ID=EGHTMLX) FILE=EGHTMLX
    OPTIONS(BITMAP_MODE='INLINE')
    %HTML5AccessibleGraphSupported
    ENCODING='utf-8'
    STYLE=HTMLBlue
    NOGTITLE
    NOGFOOTNOTE
    GPATH=&sasworklocation
;

/*   START OF NODE: spgm_db_chemo_regimen_250801.sas   */
%LET _CLIENTTASKLABEL='spgm_db_chemo_regimen_250801.sas';
%LET _CLIENTPROCESSFLOWNAME='Process Flow';
%LET _CLIENTPROJECTPATH='C:\Users\Admin\Desktop\sas_code\GangMira.egp';
%LET _CLIENTPROJECTPATHHOST='KYXU99YY10001ZF';
%LET _CLIENTPROJECTNAME='GangMira.egp';
%LET _SASPROGRAMFILE='C:\Users\Admin\Desktop\sas_code\spgm_db_chemo_regimen_250801.sas';
%LET _SASPROGRAMFILEHOST='KYXU99YY10001ZF';

/*Skim the original data*/
proc print data=CDWDATA.t100000000099034_20250801_114242(obs=5); run;
proc contents data=cdwdata.t100000000099034_20250801_114242; run;
/*Change the columns' names to comprehensible one.
Check the 'Label' printed by the 'contents' function for comparison.*/
data df_chemo;
	set cdwdata.t100000000099034_20250801_114242(
		rename=('LINE_______5'n=line '_______3'n=sex '_____________1'n=id 
						'_____________4'n=date_written '_________________10'n=diag_code 
						'_______________________6'n=regimen_code '_______________________7'n=regimen_name 
						'_____________________________8'n=date_init_regimen 
						'________________________________'n=age_init_regimen 
						'___________________cycle__9'n=cycle));
run;
/*Check the data with renamed columns*/
proc contents data=df_chemo; run;
proc print data=df_chemo(obs=5); run;
/*Create regimen_code_reduced column from regimen_code*/
data df_chemo;
	set df_chemo;
	/*Extract substring before the first '-'.*/
	regimen_code_reduced = scan(regimen_code, 1, '-');
run;
/*See the number of unique values under the regimen_code_reduced.
Each number in the range from 1 to 15 means a certain type of cancer.
E.g. 1: Gastric Cancer, 7: Breast Cancer. */
proc freq data=df_chemo;
	tables regimen_code_reduced diag_code;
run;
/*Check date-related-columns*/
proc print data=df_chemo(obs=5);
	var date_init_regimen date_written;
run;



/*Solve the ratio of patients whose diagnosis of cancer was differed more than once*/
/*1. Sort date_init_regimen inside of unique variables of id*/
proc sort data=df_chemo;
	by id date_init_regimen;
run;
/*Save the created dataset into server for validation*/
data cdwdata.df_chemo_sorted;
	set df_chemo;
run;
/*Compare the df_chemo_sorted with the original data*/
proc contents data=cdwdata.t100000000099034_20250801_114242; run;
proc contents data=cdwdata.df_chemo_sorted; run;
/*2. Create flags for checking occurrance of differing-diagnosis-event*/
data df_chemo_change_flag;
	set df_chemo;
	/*Set unique-id-groups*/
	by id;
	/*'first_diag_code', 'diag_changed_flag' columns will be retained by running this.*/
	retain first_diag_code diag_changed_flag;

	/*"first.": first element of column up to '.'*/
	/*i.e. if current id is the first id in each unique-id-group*/
	if first.id then do;
		/*Store current 'diag_code' to 'first_diag_code'*/
		first_diag_code = diag_code;
		diag_changed_flag = 0;
	end;
	/*else if current diag_code is not 'first_diag_code'*/
	else if diag_code ne first_diag_code then do;
		diag_changed_flag = 1;
	end;
	/*"last.": last element of column up to '.'*/
	/*i.e. if current id is the last element in each unique-id-group, give the output*/
	/*The output will be given for each unique-id-group because of the "by id" above.*/
	if last.id then output;
run;
/*Save the created dataset into server for validation*/
data cdwdata.df_change_flag_by_diag_code;
	set df_chemo_change_flag;
run;
/*Compare the df_change_flag_by_diag_code with the original data*/
proc contents data=cdwdata.t100000000099034_20250801_114242; run;
/*Can see the 'Observations' reduced to the number of patients*/
proc contents data=cdwdata.df_change_flag_by_diag_code; run;
/*3. Get the ratio of patients whose diagnosis of cancer was differed over all patients*/
proc sql noprint;
	select count(*) into :count_changed
	/*From the newly created data, count the number of patients with changed diagnosis*/
	from df_chemo_change_flag
	where diag_changed_flag = 1;

	/*From the df_chemo, count the number of all patients.*/
	select count(distinct id) into :n_unique_id
	from df_chemo;

	/*From the df_chemo, count the number of patients with two or more id.*/
	select count(*) into :n_multi_id
	from (
		select id
		from df_chemo
		group by id
		having count(*) >= 2	
	);
quit;
%put NOTE: Changed Diagnosis patients = &count_changed;
%put NOTE: Total patients = &n_unique_id;
%put NOTE: Patients with two or more ids = &n_multi_id;
%put NOTE: Proportion of changed Diagnosis over all = %sysevalf(&count_changed / &n_unique_id);
%put NOTE: Proportion changed Diagnosis over two or more ids = %sysevalf(&count_changed / &n_multi_id);


/*Solve the ratio of patients whose regimen(reduced) of cancer was differed more than once*/
/*1. Sort date_init_regimen inside of unique variables of id*/
proc sort data=df_chemo;
	by id date_init_regimen;
run;
/*2. Create flags for checking occurrance of differing-diagnosis-event*/
data df_change_flag_by_reg_reduced;
	set df_chemo;
	/*Set unique-id-groups*/
	by id;
	/*'first_diag_code', 'regimen_changed_flag' columns will be retained by running this.*/
	retain first_regimen_code_reduced regimen_changed_flag;

	/*"first.": first element of column up to '.'*/
	/*i.e. if current id is the first id in each unique-id-group*/
	if first.id then do;
		/*Store current 'regimen_code_reduced' to 'first_regimen_code_reduced'*/
		first_regimen_code_reduced = regimen_code_reduced;
		regimen_changed_flag = 0;
	end;
	/*else if current regimen_code_reduced is not 'first_regimen_code_reduced'*/
	else if regimen_code_reduced ne first_regimen_code_reduced then do;
		regimen_changed_flag = 1;
	end;
	/*"last.": last element of column up to '.'*/
	/*i.e. if current id is the last element in each unique-id-group, give the output*/
	/*The output will be given for each unique-id-group because of the "by id" above.*/
	if last.id then output;
run;
/*Save the created dataset into server for validation*/
data cdwdata.df_change_flag_by_reg_reduced;
	set df_change_flag_by_reg_reduced;
run;
/*Compare the df_change_flag_by_diag_code with the original data*/
proc contents data=cdwdata.t100000000099034_20250801_114242; run;
/*Can see the 'Observations' reduced to the number of patients*/
proc contents data=cdwdata.df_change_flag_by_reg_reduced; run;
/*3. Get the ratio of patients whose diagnosis of cancer was differed over all patients*/
proc sql noprint;
	select count(*) into :count_changed
	/*From the newly created data, count the number of patients with changed diagnosis*/
	from df_change_flag_by_reg_reduced
	where regimen_changed_flag = 1;

	/*From the df_chemo, count the number of patients with two or more id.*/
	select count(*) into :n_multi_id
	from (
		select id
		from df_chemo
		group by id
		having count(*) >= 2	
	);
quit;
%put NOTE: Changed Diagnosis patients = &count_changed;
%put NOTE: Total patients = &n_unique_id;
%put NOTE: Patients with two or more ids = &n_multi_id;
%put NOTE: Proportion of changed Diagnosis over all = %sysevalf(&count_changed / &n_unique_id);
%put NOTE: Proportion changed Diagnosis over two or more ids = %sysevalf(&count_changed / &n_multi_id);


/*Solve the ratio of patients whose regimen(full) of cancer was differed more than once*/
/*1. Sort date_init_regimen inside of unique variables of id*/
proc sort data=df_chemo;
	by id date_init_regimen;
run;
/*Skim the saved-sorted-data to check the right colname for the regimen code*/
proc contents data=cdwdata.df_chemo_sorted; run;
/*2. Create flags for checking occurrance of differing-diagnosis-event*/
data df_chemo_change_flag_by_reg_full;
	set df_chemo;
	/*Set unique-id-groups*/
	by id;
	/*'first_diag_code', 'regimen_changed_flag' columns will be retained by running this.*/
	retain first_regimen_code regimen_full_changed_flag;

	/*"first.": first element of column up to '.'*/
	/*i.e. if current id is the first id in each unique-id-group*/
	if first.id then do;
		/*Store current 'regimen_code' to 'first_regimen_code'*/
		first_regimen_code = regimen_code;
		regimen_full_changed_flag = 0;
	end;
	/*else if current regimen_code is not 'first_regimen_code'*/
	else if regimen_code ne first_regimen_code then do;
		regimen_full_changed_flag = 1;
	end;
	/*"last.": last element of column up to '.'*/
	/*i.e. if current id is the last element in each unique-id-group, give the output*/
	/*The output will be given for each unique-id-group because of the "by id" above.*/
	if last.id then output;
run;
/*Save the created dataset into server for validation*/
data cdwdata.df_chemo_change_flag_by_reg_full;
	set df_chemo_change_flag_by_reg_full;
run;
/*Compare the df_change_flag_by_diag_code with the original data*/
proc contents data=cdwdata.t100000000099034_20250801_114242; run;
/*Can see the 'Observations' reduced to the number of patients*/
proc contents data=cdwdata.df_chemo_change_flag_by_reg_full; run;
/*3. Get the ratio of patients whose diagnosis of cancer was differed over all patients*/
proc sql noprint;
	select count(*) into :count_changed
	/*From the newly created data, count the number of patients with changed diagnosis*/
	from df_chemo_change_flag_by_reg_full
	where regimen_full_changed_flag = 1;

	/*From the df_chemo, count the number of patients with two or more id.*/
	select count(*) into :n_multi_id
	from (
		select id
		from df_chemo
		group by id
		having count(*) >= 2	
	);
quit;
%put NOTE: Changed Diagnosis patients = &count_changed;
%put NOTE: Total patients = &n_unique_id;
%put NOTE: Patients with two or more ids = &n_multi_id;
%put NOTE: Proportion of changed Diagnosis over all = %sysevalf(&count_changed / &n_unique_id);
%put NOTE: Proportion changed Diagnosis over two or more ids = %sysevalf(&count_changed / &n_multi_id);


/*Calculating unique number of regimens by patients with 2 or more data*/
/* 1) Create a table that keeps only patients having 2 or more rows*/
proc sql;
	create table work.multiline_patients as
	select id
	from df_chemo
	group by id
	having count(*) >= 2;
quit;

/* 2) For those patients, count distinct regimen_name per patient*/
proc sql;
	create table work.patient_regimen_counts as
	/*a: an alias for df_chemo*/
	select a.id,
		/* Counts distinct non-missing regimens per id; missing counts as 0 */
		/* Using CASE inside COUNT(DISTINCT) to avoid dropping e*/
		count(distinct case when not missing(a.regimen_name) 
			then a.regimen_name end) as n_regimens
	from df_chemo as a
	/*p: an alias for multiline_patients*/
	inner join work.multiline_patients as p
		on a.id = p.id
	group by a.id;
quit;

/* 3) Build Distribution: how many patinets fall into each n_regimens bucket */
proc sql noprint;
	/* Create a table contains distribution	*/
	create table work.regimen_dist as
	select n_regimens,
		count(*) as n_patients
	from work.patient_regimen_counts
	group by n_regimens
	order by n_regimens;

	/* total number of patients in the "overall" set (denominator) */
	select count(*) into :overall_total
	from work.patient_regimen_counts;
quit;
/* Show the created distibution table */
proc print data=regimen_dist; run;
/* 4) Compute percentage and format to 0.1% */
data work.regimen_dist;
	set work.regimen_dist;
	pct = n_patients / &overall_total;
	/*
	- percentn: locale-independent format of percentage
	- 8: Maximum 8 characters are allowed to each output value
	- .2: Max number after the decimal point
	*/
	format pct percentn8.2;
run;

/* Plotting by percentage of number of different regimen */
/* --- 1) Bucket n_regimens into 1,2,3,4,>=5 and aggregate percent --- */
data work.regimen_dist;
	set work.regimen_dist;
	length bucket $6;
	if n_regimens >= 5 then bucket = '>=5';
	else bucket = strip(put(n_regimens, best.));
run;
/* Aggregate percent (and counts) by bucket */
proc sql;
	create table work.pie_data as
	select
		bucket,
		sum(n_patients) as n_patients,
		sum(pct) as pct
	from work.regimen_dist
	group by bucket
	/*Order buckets as 1, 2, 3, 4, >=5 for nicer look of plot*/
	order by case when bucket='1' then 1
						when bucket='2' then 2
						when bucket='3' then 3
						when bucket='4' then 4
						else 5 end;
quit;

/* --- 2) Pie chart --- */
title "Share by Distinct Regimen Count (1~4, >=5) among Patients with >= Rows";
proc gchart data=work.pie_data;
	pie bucket /
		sumvar=pct
		type=sum
		other=0
		midpoints=1 2 3 4 5
		percent=outside
		value=none
		slice=outside;
	format pct percentn8.2;
run;
title;

proc print data=regimen_dist
/* --- 3) Barplot --- */
title "Percent by Distinct Regimen Count";
proc sgplot data=work.regimen_dist;
	vbarparm category=n_regimens response=pct / datalabel;
	xaxis label="Distinct regimen count (n_regimens)" discreteorder=data;
	yaxis label="Percent of overall" valuesformat=percentn8.2;
	format pct percentn8.2;
run;
title;



/* Count number of discrepancy in `regimen_name` by `line`  */
/* Count number of data grouped by the 'id' column	*/
/* In 'PROC SUMMARY' or 'PROC MEANS', 
numbers of rows grouped by class will be automatically goes into _FREQ_ variable. */
proc summary data=cdwdata.df_chemo_sorted nway;
	class id;
	output out=_count_by_pt(drop=_type_) / autoname;
run;
/*1)  Create a data with 2 or more ids per patient */
data multi_ids;
	set _count_by_pt;
	if _freq_ >= 2;
	keep id;
run;
proc contents data=multi_ids; run;
/* Filter data from df_chemo with multi_ids */
proc sort data=multi_ids; by id; run;
/*Change the columns' names to comprehensible one.
Check the 'Label' printed by the 'contents' function for comparison.*/
data df_chemo;
	set cdwdata.t100000000099034_20250801_114242(
		rename=('LINE_______5'n=line '_______3'n=sex '_____________1'n=id 
						'_____________4'n=date_written '_________________10'n=diag_code 
						'_______________________6'n=regimen_code '_______________________7'n=regimen_name 
						'_____________________________8'n=date_init_regimen 
						'________________________________'n=age_init_regimen 
						'___________________cycle__9'n=cycle));
run;
proc sort data=df_chemo;
	by id date_init_regimen;
run;
/* Obs = 134169 : number of the entire data */
proc contents data=df_chemo; run; 
data df_multi;
	merge df_chemo(in=a) multi_ids(in=b);
	by id;
	if a and b;
run;
/* Obs = 103256 : number of data only with multi-cohort*/
proc contents data=df_multi; run;
/* 2) Drop NAs in regimen_name */
data df_multi_no_na;
	set df_multi;
	if not missing(regimen_name);
run;
/* Obs = 103160 : number of data only with multi-cohort after dropping NAs */
proc contents data=df_multi_no_na; run;
/* 3) Make unique triplets (id, line, regimen_name) by sorting df_multi_no_na */
proc sort data=df_multi_no_na out=_dedup_reg nodupkey;
	by id line regimen_name;
run;
/* 4) Count unique regimen per (id, line) */
proc freq data=_dedup_reg noprint;
	tables id*line / out=_per_line_unique(drop=percent);
run;
/* rename COUNT -> n_reg_unique in '_dedup_reg' */
data _per_line_unique;
	set _per_line_unique;
	n_reg_unique = count;
	drop count;
run;
/* Obs = 89074 : number of data after filtering unique regimens per patient */
proc contents data=_per_line_unique; run;
proc print data=_per_line_unique(obs=30); run;
/* 5) discrepancy per line := n_reg_unique - 1 */
/* Compare the outcome with the DF_CHEMO_SORTED for the clarification */
data _per_line_discrep;
	set _per_line_unique;
	/* Use max() for prevent any potential negative value */
	discrep = max(0, n_reg_unique - 1);
run;
/* 6) Sum by patient -> df_line_regimen_discrep (fill NA with 0) */
proc summary data=_per_line_discrep nway;
	class id;
	var discrep;
	output out=_sum_by_pt(drop=_type_ _freq_) sum(discrep)=n_discrep;
run;
/* Sort the _sum_by_pt by id to merge with multi_ids data*/
proc sort data=_sum_by_pt; 
	by id; 
run;
data df_line_regimen_discrep;
	merge multi_ids(in=b) _sum_by_pt(in=c);
	by id;
	/* Keep multi cohort only */
	if b;
	/* Patients with no discrepancy won't be in the _sum_by_pt */
	if not c then n_discrep = 0;
run;
/* 7) data preparation for a pie chart */
data _null_;
	if 0 then set df_line_regimen_discrep nobs=n;
	call symputx('TOTAL_MULTI', n);
	stop;
run;
proc summary data=df_line_regimen_discrep nway;
	class n_discrep;
	output out=discrep_freq(drop=_type_) / autoname;
run;
data discrep_freq;
	set discrep_freq;
	pct = _freq_ / &total_multi.;
	format pct percentn8.2;
	length label $12;
	/* Create a character column 'label' for the chart */
	label = put(pct, percentn8.2);
run;
proc contents data=discrep_freq; run;
/* obs = 32949 : the number of all patients with multi-cohort */
proc print data=discrep_freq; var; sum _freq_; run;
/* Check if the upper number is correct */
proc contents data=multi_ids; run;
/* Count the number of patients with at least 1 discrepancy: 11551*/
proc sql;
	select sum(_freq_)
	from discrep_freq
	where n_discrep >= 1;
quit;
data work.discrep_freq;
	set discrep_freq;
run;
/* 8) Draw the pie chart*/
/* --- 1) Bucket n_discrep into 1,2,3,4>= and aggregate percent --- */
data work.discrep_freq_dist;
	set discrep_freq;
	/* $3: 3 characters are enough for representing '4>=' */
	length bucket $3;
	if n_discrep >= 4 then bucket = '>=4';
	else bucket = strip(put(n_discrep, best.));
run;
/* Aggregate percent (and counts) by bucket */
proc sql;
	create table work.pie_discrep_freq as
	select
		bucket,
		sum(_freq_) as n_patients,
		sum(pct) as pct
	from work.discrep_freq_dist
	group by bucket
	/*Order buckets as 1, 2, 3, >=4 for nicer look of plot*/
	order by case when bucket='0' then 1
						when bucket='1' then 2
						when bucket='2' then 3
						when bucket='3' then 4
						else 5 end;
quit;
title "Discrepancies per Patient with Multi-regimen Cohort";
proc gchart data=pie_discrep_freq;
	pie bucket / 
		sumvar=pct
		type=sum
		other=0
		midpoints=1 2 3 4 5
		percent=outside
		value=none
		slice=outside;
	format pct percentn8.2;
run;
title;

title "Percent by number of Discrepancy";
proc sgplot data=discrep_freq;
	vbarparm category=n_discrep response=pct / datalabel;
	xaxis label="Discrepancy count (n_discrep)" discreteorder=data;
	yaxis label="Percent of overall" valuesformat=percentn8.2;
	format pct percentn8.2;
run;
title;

%LET _CLIENTTASKLABEL=;
%LET _CLIENTPROCESSFLOWNAME=;
%LET _CLIENTPROJECTPATH=;
%LET _CLIENTPROJECTPATHHOST=;
%LET _CLIENTPROJECTNAME=;
%LET _SASPROGRAMFILE=;
%LET _SASPROGRAMFILEHOST=;

;*';*";*/;quit;run;
ODS _ALL_ CLOSE;
