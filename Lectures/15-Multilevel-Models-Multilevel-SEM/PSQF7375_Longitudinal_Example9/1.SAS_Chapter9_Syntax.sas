* Output options: all can be turned on/off by adding or removing the NO;
* page number, date, centering, or page breaks, page length and width;
OPTIONS NOnumber NOdate NOcenter FormDlim=' ' PageSize=MAX LineSize=MAX;
* Eliminate SAS default titles and names of tables in output (TRACE ON to show);
TITLE; ODS TRACE OFF; 

* Defining global variable for file location to be replaced in code below;
%LET filesave=C:\Dropbox\19_PSQF7375_Longitudinal\PSQF7375_Longitudinal_Example9;
* Location for SAS files for these models (uses macro variable filesave);
LIBNAME filesave "&filesave.";

* Import chapter 9 multivariate data in work library;
* Also make person mean monitoring;
DATA work.Chapter9multiv; SET filesave.SAS_Chapter9;
PMmonitor = MEAN(OF monitor12-monitor18);
LABEL PMmonitor = "PMmonitor: Person Mean Monitoring";
RUN;

* Stack multivariate data;
* Also make a copy of monitor18 for later;
DATA work.Chapter9; SET work.Chapter9multiv;
copymonitor18 = monitor18;
occasion=12; age=age12; risky=risky12; monitor=monitor12; OUTPUT;
occasion=13; age=age13; risky=risky13; monitor=monitor13; OUTPUT;
occasion=14; age=age14; risky=risky14; monitor=monitor14; OUTPUT;
occasion=15; age=age15; risky=risky15; monitor=monitor15; OUTPUT;
occasion=16; age=age16; risky=risky16; monitor=monitor16; OUTPUT;
occasion=17; age=age17; risky=risky17; monitor=monitor17; OUTPUT;
occasion=18; age=age18; risky=risky18; monitor=monitor18; OUTPUT;
* Drop old unnecessary multivariate variables;
DROP age12-age18 risky12-risky18 monitor12-monitor18;
LABEL   
copymonitor18 = "copymonitor18: Monitoring at Age 18 (copy)"
occasion = "occasion: Occasion of Measurement (12-18)"
age = "age: Exact Age at Occasion"
risky = "risky: Risky Behavior at Occasion"
monitor = "monitor: Monitoring at Occasion";
RUN;

* Center predictors for analysis in stacked data;
DATA work.Chapter9; SET work.Chapter9;
agec18 = age - 18;
att4 = attitude12 - 4;
mon3 = monitor - 3;
PMmon3 = PMmonitor - 3;
WPmon = monitor - PMmonitor;
Age18mon3 = copymonitor18 - 3;
Change18mon = monitor - copymonitor18;
LABEL
agec18 = "agec18: Exact Age (0=18)"
att4 = "att4: Age 12 Attitudes (0=4)"
mon3 = "mon3: Monitoring (0=3)"
PMmon3 = "PMmon3: Person Mean Monitoring (0=3)"
WPmon = "WPmon: Within-Person Monitoring (0=PM)"
Age18mon3 = "Age18mon3: BP Monitoring at Age 18 (0=3)"
Change18mon = "Change18mon: WP Monitoring from Age 18 (0=Age18)";
* Subset sample to complete cases for all predictors;
IF NMISS(agec18, att4, risky, PMmonitor, Age18mon3, monitor)>0 THEN DELETE;
RUN;
 
* Trimming data to send just needed variables to Mplus;
DATA work.Chapter9Mplus; SET work.Chapter9; 
agesq=agec18*agec18;
* Telling it which variables to keep -- handy to use this in Mplus;
KEEP PersonID occasion risky agec18 att4 mon3 agesq;
* Replace any missing values with -999;
	ARRAY avars(7) FamilyID occasion risky agec18 att4 agesq mon3;
		DO i=1 TO 7; IF avars(i)=. THEN avars(i)=-999; END; DROP i;
RUN;

* Export to .csv for use in Mplus MLM syntax;
PROC EXPORT DATA=work.Chapter9Mplus OUTFILE= "&filesave.\Chapter9.csv" 
	DBMS=CSV REPLACE; PUTNAMES=NO; RUN;


* Create double-stacked dataset for multivariate analysis;
* Also create dummy codes to serve as intercepts if needed;
DATA work.Chapter9doublestack; SET work.Chapter9;
DV="1risky  "; dvR=1; dvM=0; outcome=risky; OUTPUT;
DV="2monitor"; dvR=0; dvM=1; outcome=mon3;  OUTPUT;
LABEL
DV = "DV: Categorical indicator for which DV the row is for"
dvR = "dvR: Intercept 1=risky, 0=monitor"
dvM = "dvM: Intercept 0=risky, 1=monitor"
outcome = "outcome: Combined outcome variable column";
RUN;

* Open output directory to save results to;
ODS RTF FILE="&filesave.\Example9_SAS_Output.rtf" STYLE=HTMLBlue;

TITLE1 'Ch 9 Eq 9.3: Multivariate Model of Risky Behavior and Monitoring';
TITLE2 'Tricking Univariate Software';
PROC MIXED DATA=work.Chapter9doublestack COVTEST NOCLPRINT NAMELEN=100 IC METHOD=ML;
     CLASS PersonID occasion DV;
     MODEL outcome = dvR dvM dvR*agec18 dvM*agec18 dvR*agec18*agec18 
                     dvR*att4 dvR*agec18*att4
                       / NOINT SOLUTION DDFM=Satterthwaite;
     RANDOM dvR dvM dvR*agec18 dvM*agec18 / G GCORR TYPE=UN SUBJECT=PersonID;
     REPEATED DV / R RCORR TYPE=UN SUBJECT=PersonID*occasion;
RUN; TITLE1; TITLE2;

* Close output directory;
ODS RTF CLOSE;

