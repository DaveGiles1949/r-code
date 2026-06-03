# LOADING OF DATA BY YEAR,AND BASIC PLOTS
# INCLUDING PLOTS OF THE BENFORD ANALYSIS OF THE 1ST. AND 2ND. DIGITS OF THE COEFFICIENTS,
# STD. ERRORS, AND T-STATISTICS
# =======================================

library(benford.analysis)
library(stringr)
library(readxl)
c<- c()
s<- c()
t<- c()

coef<- read_excel("C:/Users/OEM/Sync/Benford Law/NZEP/NZEP.xlsx",sheet="coef", col_names=TRUE)
coef_1966<- coef$coef66[!is.na(coef$coef66)]
coef_1968<- coef$coef68[!is.na(coef$coef68)]
coef_1969<- coef$coef69[!is.na(coef$coef69)]
coef_1970<- coef$coef70[!is.na(coef$coef70)]
coef_1971<- coef$coef71[!is.na(coef$coef71)]
coef_1972<- coef$coef72[!is.na(coef$coef72)]
coef_1973<- coef$coef73[!is.na(coef$coef73)]
coef_1974<- coef$coef74[!is.na(coef$coef74)]
coef_1975<- coef$coef75[!is.na(coef$coef75)]
coef_1976<- coef$coef76[!is.na(coef$coef76)]
coef_1977<- coef$coef77[!is.na(coef$coef77)]
coef_1978<- coef$coef78[!is.na(coef$coef78)]
coef_1979<- coef$coef79[!is.na(coef$coef79)]
coef_1980<- coef$coef80[!is.na(coef$coef80)]
coef_1981<- coef$coef81[!is.na(coef$coef81)]
# No regressions published in 1982
coef_1983<- coef$coef83[!is.na(coef$coef83)]
# No regressions published in 1984
coef_1985<- coef$coef85[!is.na(coef$coef85)]
coef_1986<- coef$coef86[!is.na(coef$coef86)]
coef_1987<- coef$coef87[!is.na(coef$coef87)]
coef_1988<- coef$coef88[!is.na(coef$coef88)]
coef_1989<- coef$coef89[!is.na(coef$coef89)]
coef_1990<- coef$coef90[!is.na(coef$coef90)]
coef_1991<- coef$coef91[!is.na(coef$coef91)]
coef_1992<- coef$coef92[!is.na(coef$coef92)]
coef_1993<- coef$coef93[!is.na(coef$coef93)]
coef_1994<- coef$coef94[!is.na(coef$coef94)]
coef_1995<- coef$coef95[!is.na(coef$coef95)]
coef_1996<- coef$coef96[!is.na(coef$coef96)]
coef_1997<- coef$coef97[!is.na(coef$coef97)]
coef_1998<- coef$coef98[!is.na(coef$coef98)]
coef_1999<- coef$coef99[!is.na(coef$coef99)]
coef_2000<- coef$coef00[!is.na(coef$coef00)]
coef_2001<- coef$coef01[!is.na(coef$coef01)]
coef_2002<- coef$coef02[!is.na(coef$coef02)]
coef_2003<- coef$coef03[!is.na(coef$coef03)]
coef_2004<- coef$coef04[!is.na(coef$coef04)]
coef_2005<- coef$coef05[!is.na(coef$coef05)]
coef_2006<- coef$coef06[!is.na(coef$coef06)]
coef_2007<- coef$coef07[!is.na(coef$coef07)]
coef_2008<- coef$coef08[!is.na(coef$coef08)]
coef_2009<- coef$coef09[!is.na(coef$coef09)]
coef_2010<- coef$coef10[!is.na(coef$coef10)]
coef_2011<- coef$coef11[!is.na(coef$coef11)]
coef_2012<- coef$coef12[!is.na(coef$coef12)]
coef_2013<- coef$coef13[!is.na(coef$coef13)]
coef_2014<- coef$coef14[!is.na(coef$coef14)]
coef_2015<- coef$coef15[!is.na(coef$coef15)]
coef_2016<- coef$coef16[!is.na(coef$coef16)]
coef_2017<- coef$coef17[!is.na(coef$coef17)]
coef_2018<- coef$coef18[!is.na(coef$coef18)]
coef_2019<- coef$coef19[!is.na(coef$coef19)]
coef_2020<- coef$coef20[!is.na(coef$coef20)]
coef_2021<- coef$coef21[!is.na(coef$coef21)]
coef_2022<- coef$coef22[!is.na(coef$coef22)]
coef_2023<- coef$coef23[!is.na(coef$coef23)]
coef_2024<- coef$coef24[!is.na(coef$coef24)]
coef_2025<- coef$coef25[!is.na(coef$coef25)]

coef_total<- c(coef_1966, coef_1968,coef_1969,coef_1970,coef_1971,coef_1972,coef_1973,coef_1974,
coef_1975,coef_1976,coef_1977,coef_1978,coef_1979,coef_1980,coef_1981,coef_1983,coef_1985,coef_1986,
coef_1987,coef_1988,coef_1989,coef_1990,coef_1991,coef_1992,coef_1993,coef_1994,coef_1995,
coef_1996,coef_1997,coef_1998,coef_1999,coef_2000,coef_2001,coef_2002,coef_2003,coef_2004,coef_2005,
coef_2006,coef_2007,coef_2008,coef_2009,coef_2010,coef_2011,coef_2012,coef_2013,coef_2014,coef_2015,coef_2016,
coef_2017,coef_2018,coef_2019,coef_2020,coef_2021,coef_2022,coef_2023,coef_2024,coef_2025)

# NOTE: The ordering of the 2024 std errors differs from the ordering of their associated coefficients (within the year)
# but this doesn't affect any calculations

se<- read_excel("C:/Users/OEM/Sync/Benford Law/NZEP/NZEP.xlsx",sheet="se", col_names=TRUE)
se_1966<- se$se66[!is.na(se$se66)]
se_1968<- se$se68[!is.na(se$se68)]
se_1969<- se$se69[!is.na(se$se69)]
# No s.e.'s in 1970
se_1971<- se$se71[!is.na(se$se71)]
se_1972<- se$se72[!is.na(se$se72)]
se_1973<- se$se73[!is.na(se$se73)]
se_1974<- se$se74[!is.na(se$se74)]
se_1975<- se$se75[!is.na(se$se75)]
# No s.e.'s in 1976
se_1977<- se$se77[!is.na(se$se77)]
# No s.e.'s in 1978
se_1979<- se$se79[!is.na(se$se79)]
se_1980<- se$se80[!is.na(se$se80)]
se_1981<- se$se81[!is.na(se$se81)]
# No regressions published in 1982
se_1983<- se$se83[!is.na(se$se83)]
# No regressions published in 1984
# No s.e.'s in 1985-1987
se_1988<- se$se88[!is.na(se$se88)]
# No s.e.'s in 1989
se_1990<- se$se90[!is.na(se$se90)]
# No s.e.'s in 1991
se_1992<- se$se92[!is.na(se$se92)]
# No s.e.'s in 1993 or 1994
se_1995<- se$se95[!is.na(se$se95)]
se_1996<- se$se96[!is.na(se$se96)]
se_1997<- se$se97[!is.na(se$se97)]
se_1998<- se$se98[!is.na(se$se98)]
se_1999<- se$se99[!is.na(se$se99)]
se_2000<- se$se00[!is.na(se$se00)]
se_2001<- se$se01[!is.na(se$se01)]
# No s.e.'s in 2002
se_2003<- se$se03[!is.na(se$se03)]
se_2004<- se$se04[!is.na(se$se04)]
se_2005<- se$se05[!is.na(se$se05)]
se_2006<- se$se06[!is.na(se$se06)]
se_2007<- se$se07[!is.na(se$se07)]
se_2008<- se$se08[!is.na(se$se08)]
se_2009<- se$se09[!is.na(se$se09)]
se_2010<- se$se10[!is.na(se$se10)]
se_2011<- se$se11[!is.na(se$se11)]
se_2012<- se$se12[!is.na(se$se12)]
se_2013<- se$se13[!is.na(se$se13)]
se_2014<- se$se14[!is.na(se$se14)]
se_2015<- se$se15[!is.na(se$se15)]
se_2016<- se$se16[!is.na(se$se16)]
se_2017<- se$se17[!is.na(se$se17)]
se_2018<- se$se18[!is.na(se$se18)]
se_2019<- se$se19[!is.na(se$se19)]
se_2020<- se$se20[!is.na(se$se20)]
se_2021<- se$se21[!is.na(se$se21)]
se_2022<- se$se22[!is.na(se$se22)]
se_2023<- se$se23[!is.na(se$se23)]
se_2024<- se$se24[!is.na(se$se24)]
se_2025<- se$se25[!is.na(se$se25)]

se_total<- c(se_1966,se_1968,se_1969,se_1971,se_1972,se_1973,se_1974,se_1975,se_1977,se_1979,
se_1980,se_1981,se_1983,se_1988,se_1990,se_1992,se_1995,se_1996,se_1997,se_1998,
se_1999,se_2000,se_2001,se_2003,se_2004,se_2005,se_2006,se_2007,se_2008,se_2009,se_2010,
se_2011,se_2012,se_2013,se_2014,se_2015,se_2016,se_2017,se_2018,se_2019,se_2020,se_2021,se_2022,se_2023,se_2024,se_2025)

# Read the t-stats:

tstat<- read_excel("C:/Users/OEM/Sync/Benford Law/NZEP/NZEP.xlsx",sheet="tstat", col_names=TRUE)
# no t-stats in 1966 , 1968, or 1969
tstat_1970<- tstat$tstat70[!is.na(tstat$tstat70)]
tstat_1971<- tstat$tstat71[!is.na(tstat$tstat71)]
# no t-stats in 1972
tstat_1973<- tstat$tstat73[!is.na(tstat$tstat73)]
tstat_1974<- tstat$tstat74[!is.na(tstat$tstat74)]
tstat_1975<- tstat$tstat75[!is.na(tstat$tstat75)]
tstat_1976<- tstat$tstat76[!is.na(tstat$tstat76)]
tstat_1977<- tstat$tstat77[!is.na(tstat$tstat77)]
tstat_1978<- tstat$tstat78[!is.na(tstat$tstat78)]
tstat_1979<- tstat$tstat79[!is.na(tstat$tstat79)]
# no t-stats in 1980 or 1981; no regression results in 1982
tstat_1983<- tstat$tstat83[!is.na(tstat$tstat83)]
# No regression results in 1984
tstat_1985<- tstat$tstat85[!is.na(tstat$tstat85)]
tstat_1986<- tstat$tstat86[!is.na(tstat$tstat86)]
tstat_1987<- tstat$tstat87[!is.na(tstat$tstat87)]
tstat_1988<- tstat$tstat88[!is.na(tstat$tstat88)]
tstat_1989<- tstat$tstat89[!is.na(tstat$tstat89)]
tstat_1990<- tstat$tstat90[!is.na(tstat$tstat90)]
tstat_1991<- tstat$tstat91[!is.na(tstat$tstat91)]
tstat_1992<- tstat$tstat92[!is.na(tstat$tstat92)]
tstat_1993<- tstat$tstat93[!is.na(tstat$tstat93)]
tstat_1994<- tstat$tstat94[!is.na(tstat$tstat94)]
tstat_1995<- tstat$tstat95[!is.na(tstat$tstat95)]
tstat_1996<- tstat$tstat96[!is.na(tstat$tstat96)]
tstat_1997<- tstat$tstat97[!is.na(tstat$tstat97)]
tstat_1998<- tstat$tstat98[!is.na(tstat$tstat98)]
tstat_1999<- tstat$tstat99[!is.na(tstat$tstat99)]
tstat_2000<- tstat$tstat00[!is.na(tstat$tstat00)]
tstat_2001<- tstat$tstat01[!is.na(tstat$tstat01)]
tstat_2002<- tstat$tstat02[!is.na(tstat$tstat02)]
tstat_2003<- tstat$tstat03[!is.na(tstat$tstat03)]
tstat_2004<- tstat$tstat04[!is.na(tstat$tstat04)]
# no t-stats in 2005
tstat_2006<- tstat$tstat06[!is.na(tstat$tstat06)]
# no t-stats in 2007
tstat_2008<- tstat$tstat08[!is.na(tstat$tstat08)]
# no t-stats in 2009
tstat_2010<- tstat$tstat10[!is.na(tstat$tstat10)]
tstat_2011<- tstat$tstat11[!is.na(tstat$tstat11)]
tstat_2012<- tstat$tstat12[!is.na(tstat$tstat12)]
# no t-stats in 2013
tstat_2014<- tstat$tstat14[!is.na(tstat$tstat14)]
tstat_2015<- tstat$tstat15[!is.na(tstat$tstat15)]
tstat_2016<- tstat$tstat16[!is.na(tstat$tstat16)]
tstat_2017<- tstat$tstat17[!is.na(tstat$tstat17)]
tstat_2018<- tstat$tstat18[!is.na(tstat$tstat18)]
tstat_2019<- tstat$tstat19[!is.na(tstat$tstat19)]
# no t-stats in 2020
tstat_2021<- tstat$tstat21[!is.na(tstat$tstat21)]
tstat_2022<- tstat$tstat22[!is.na(tstat$tstat22)]
# no t-stats in 2023
tstat_2024<- tstat$tstat24[!is.na(tstat$tstat24)]
# no t-stats in 2025

tstat_total<- c(tstat_1970,tstat_1971, tstat_1973,tstat_1974,tstat_1975,tstat_1976,tstat_1977,tstat_1978,tstat_1979,
tstat_1983,tstat_1985,tstat_1986,tstat_1987,tstat_1988,tstat_1989,tstat_1990,tstat_1991,tstat_1992,
tstat_1993,tstat_1994,tstat_1995,tstat_1996,tstat_1997,tstat_1998,tstat_1999,tstat_2000,
tstat_2001,tstat_2002,tstat_2003,tstat_2004,tstat_2006,tstat_2008,tstat_2010,tstat_2011,tstat_2012,tstat_2014,
tstat_2015,tstat_2016,tstat_2017,tstat_2018,tstat_2019,tstat_2021,tstat_2022,tstat_2024)

# In most years, the t-stat. data on file are absolute values. In some years the negative signs were recorded
# No attention should be paid to signs associated with t-statistics in the data
# Absolute values of the t-stats are created when extracting the first 2 digits below, 
# So, all t-stat. values are positive in the actual Benford analysis

vol<- seq(1, 59,1)
year<- seq(1966,2025,1)
# NOTE: Volume 1 of "NZEP" was in 1966; Volume 2 was in 1968. There was no volume in 1967.

# Np = total no. of papers published per volume
# np = no. of papers per volume with regression results analyzed in this study

Np<- c(10,10,11,7,10,7,8,8,11,9,10,9,12,8,10,9,15,9,9,10,14,5,4,10,14,10,12,
12,12,12,12,11,13,13,11,25,12,14,8,8,11,8,16,14,18,20,17,23,16,16,19,19,17, 15,19,27,28,18,21)

np<- c(3,2,1,1,5,2,5,4,7,3,2,3,3,5,2,0,5,0,1,3,2,2,1,3,2,5,3,4,5,5,5,6,7,4,4,1,3,5,2,3,6,4, 5,
5,9,9,7,8, 8, 6,6,13,8,10,12,11,4,14,6)

summary(np/Np)
sd(np/Np)

c[59]<- length(coef_2025)
c[58]<- length(coef_2024)
c[57]<- length(coef_2023)
c[56]<- length(coef_2022)
c[55]<- length(coef_2021)
c[54]<- length(coef_2020)
c[53]<- length(coef_2019)
c[52]<- length(coef_2018)
c[51]<- length(coef_2017)
c[50]<- length(coef_2016)
c[49]<- length(coef_2015)
c[48]<- length(coef_2014)
c[47]<- length(coef_2013)
c[46]<- length(coef_2012)
c[45]<- length(coef_2011)
c[44]<- length(coef_2010)
c[43]<- length(coef_2009)
c[42]<- length(coef_2008)
c[41]<- length(coef_2007)
c[40]<- length(coef_2006)
c[39]<- length(coef_2005)
c[38]<- length(coef_2004)
c[37]<- length(coef_2003)
c[36]<- length(coef_2002)
c[35]<- length(coef_2001)
c[34]<- length(coef_2000)
c[33]<- length(coef_1999)
c[32]<- length(coef_1998)
c[31]<- length(coef_1997)
c[30]<- length(coef_1996)
c[29]<- length(coef_1995)
c[28]<- length(coef_1994)
c[27]<- length(coef_1993)
c[26]<- length(coef_1992)
c[25]<- length(coef_1991)
c[24]<- length(coef_1990)
c[23]<- length(coef_1989)
c[22]<- length(coef_1988)
c[21]<- length(coef_1987)
c[20]<- length(coef_1986)
c[19]<- length(coef_1985)
c[18]<- 0                 # There were no papers within regression rssults in 1982 or 1984
c[17]<- length(coef_1983)
c[16]<- 0
c[15]<- length(coef_1981)
c[14]<- length(coef_1980)
c[13]<- length(coef_1979)
c[12]<- length(coef_1978)
c[11]<- length(coef_1977)
c[10]<- length(coef_1976)
c[9]<- length(coef_1975)
c[8]<- length(coef_1974)
c[7]<- length(coef_1973)
c[6]<- length(coef_1972)
c[5]<- length(coef_1971)
c[4]<- length(coef_1970)
c[3]<- length(coef_1969)
c[2]<- length(coef_1968)
c[1]<- length(coef_1966)
					# In some cases, coeffcients are reported without either a std. error or a t-stat.
					# Instead, a p-value might be reported, of the level of signifiance may be flagged (*,**, etc.)
					# p-values have not been analyzed
					# If a paper reports BOTH a std. error and t-statistic for the SAME coefficient, just the s.e. is used
					# In some cases a paper reports s.e.'s for some coeffs., and t-stats. for other coeffs.		
					# In this case, the SEPARATE s.e.'s and t-stats. are used with their respective coeffs.
s[59]<- length(se_2025)
s[58]<- length(se_2024)
s[57]<- length(se_2023)
s[56]<- length(se_2022)
s[55]<- length(se_2021)
s[54]<- length(se_2020)
s[53]<- length(se_2019)
s[52]<- length(se_2018)
s[51]<- length(se_2017)
s[50]<- length(se_2016)
s[49]<- length(se_2015)
s[48]<- length(se_2014)
s[47]<- length(se_2013)
s[46]<- length(se_2012)
s[45]<- length(se_2011)
s[44]<- length(se_2010)
s[43]<- length(se_2009)
s[42]<- length(se_2008)
s[41]<- length(se_2007)
s[40]<- length(se_2006)
s[39]<- length(se_2005)
s[38]<- length(se_2004)
s[37]<- length(se_2003)
s[36]<- 0
s[35]<- length(se_2001)
s[34]<- length(se_2000)
s[33]<- length(se_1999)
s[32]<- length(se_1998)
s[31]<- length(se_1997)
s[30]<- length(se_1996)
s[29]<- length(se_1995)
s[28]<- 0
s[27]<- 0
s[26]<- length(se_1992)
s[25]<- 0
s[24]<- length(se_1990)
s[23]<- 0
s[22]<- length(se_1988)
s[21]<- 0
s[20]<- 0
s[19]<- 0
s[18]<- 0
s[17]<- length(se_1983)
s[16]<- 0
s[15]<- length(se_1981)
s[14]<- length(se_1980)
s[13]<- length(se_1979)
s[12]<- 0
s[11]<- length(se_1977)
s[10]<- 0
s[9]<- length(se_1975)
s[8]<- length(se_1974)
s[7]<- length(se_1973)
s[6]<- length(se_1972)
s[5]<- length(se_1971)
s[4]<- 0
s[3]<- length(se_1969)
s[2]<- length(se_1968)
s[1]<- length(se_1966)
      				# t-statistics are used only if they were reported as such. 
                             	#They are not computed from the coeff. and s.e. values becuase of rounding effectson digits

t[59]<- 0  
t[58]<- length(tstat_2024)
t[57]<- 0
t[56]<- length(tstat_2022)
t[55]<- length(tstat_2021)
t[54]<- 0
t[53]<- length(tstat_2019)
t[52]<- length(tstat_2018)
t[51]<- length(tstat_2017)
t[50]<- length(tstat_2016)
t[49]<- length(tstat_2015)
t[48]<- length(tstat_2014)
t[47]<- 0
t[46]<- length(tstat_2012)
t[45]<- length(tstat_2011)
t[44]<- length(tstat_2010)
t[43]<- 0
t[42]<- length(tstat_2008)
t[41]<- 0
t[40]<- length(tstat_2006)
t[39]<- 0
t[38]<- length(tstat_2004)
t[37]<- length(tstat_2003)
t[36]<- length(tstat_2002)
t[35]<- length(tstat_2001)
t[34]<- length(tstat_2000)
t[33]<- length(tstat_1999)
t[32]<- length(tstat_1998)
t[31]<- length(tstat_1997)
t[30]<- length(tstat_1996)
t[29]<- length(tstat_1995)
t[28]<- length(tstat_1994)
t[27]<- length(tstat_1993)
t[26]<- length(tstat_1992)
t[25]<- length(tstat_1991)
t[24]<- length(tstat_1990)
t[23]<- length(tstat_1989)
t[22]<- length(tstat_1988)
t[21]<- length(tstat_1987)
t[20]<- length(tstat_1986)
t[19]<- length(tstat_1985)
t[18]<- 0
t[17]<- length(tstat_1983)
t[16]<- 0
t[15]<- 0
t[14]<- 0
t[13]<- length(tstat_1979)
t[12]<- length(tstat_1978)
t[11]<- length(tstat_1977)
t[10]<- length(tstat_1976)
t[9]<- length(tstat_1975)
t[8]<- length(tstat_1974)
t[7]<- length(tstat_1973)
t[6]<- 0
t[5]<- length(tstat_1971)
t[4]<- length(tstat_1970)
t[3]<- 0 
t[2]<- 0
t[1]<- 0

summary(coef_total)
summary(se_total)
summary(tstat_total)
summary(c)
summary(s)
summary(t)
sd(c)
sd(s)
sd(t)

summary(c/np)
summary(s/np)
summary(t/np)
# Remove the "NAs"
cnp <- c/np
snp<- s/np
tnp<- t/np
cnp<- cnp[!is.na(cnp)]
snp<- snp[!is.na(snp)]
tnp<- tnp[!is.na(tnp)]
sd(cnp)
sd(snp)
sd(tnp)


# Start plotting:
################

# A plot for Benford's distributions of 1st to 3rd. digits
# ---------------------------------------------------------
par(mfrow=c(1,1))
d1<- 1:9
d2<- d3<- 0:9
zero<- rep(0,100)
benford_2<- c(0.1197,0.1139,0.1088,0.1043,0.1,0.0967, 0.0934,0.0904,0.0876,0.0850)
benford_3<- c(0.1018,0.1014,0.1010,0.1006,0.1002,0.0998,0.0994,0.0990,0.0986,0.0983)

bfdc<- benford(coef_total,1,sign="both")
bfds<- benford(se_total,1, sign="both")
bfdt<- benford(tstat_total,1, sign="both")

x1 <- seq(1,9)
x2<-0:9
y1_top <- bfdc$bfd$benford.dist
y2_top<- benford_2
y3_top<- benford_3
y1_bottom <- y2_bottom<- y3_bottom<- zero

# Create empty plot
plot(0, 0, xlim=c(0, 9), ylim=c(0, 0.35), xaxp=c(0,9,9),xlab="Digits", ylab="Relative frequency",type="n", 
cex.main=0.9,cex.lab=0.9,cex.axis=0.8,main="Figure 1:
Benford's Distributions for 1st. to 3rd. Digits")
# Add vertical line
segments(x0=x1, y0=y1_bottom, x1=x1, y1=y1_top, col="black", lwd=2)
# Add point on top
points(x1, y1_top, pch=19, col="red", cex=1.2)
#segments(x0=x2, y0=y2_bottom, x1=x2, y1=y2_top, col="black", lwd=2)
points(x2, y2_top, pch=8, col="blue", cex=1)
#segments(x0=x2, y0=y3_bottom, x1=x2, y1=y3_top, col="black", lwd=2)
points(x2, y3_top, pch=17, col="black", cex=1.2)
legend(6.5,0.3, legend=c("1st. digit", "2nd. digit", "3rd. digit"), pch=c(19,8,17),col=c("red","blue","black"), cex=0.9)

# Number of analyzed papers per issue of "NZEP"; and No. of digits analyzed per analyzed paper 
# --------------------------------------------------------------------------------------------
#par(mfrow=c(2,1))
#plot(vol,np, xlab="NZEP Volume No.", ylab="No. of Papers",, ylim=c(0,15),
#main="Figure 2:
#Number of Papers With Data per Volume", type="o",cex.main=0.8,cex.lab=0.9,cex.axis=0.8)

plot(vol,c, xlab="NZEP Volume No.", main="Figure 3:
Number of Coefficients, Standard Errors,
& t-Statistics per Volume",ylab="Sample Size", col="black", 
type="l", lty=1,lwd=1,ylim=c(1,2400),xlim=c(1,59),xaxp = c(1, 61, 30),yaxp=c(0,2400,6),cex.main=0.9,cex.lab=0.9,cex.axis=0.8)
lines(vol,s,col="red",type="l", lty=2,lwd=1)
lines(vol,t, col="blue",type="l", lty=3,lwd=1)
legend(1,2150,legend=c("coefficients", "standard errors", "t-statistics"), col=c("black","red","blue"),
lwd=c(1,1,1), lty=c(1,2,3), cex=0.8)

# Main Benford analysis for 1st. and 2nd. digits
# ----------------------------------------------
#par(mfrow=c(2,1))
plot(d1,bfdc$bfd$data.dist, col="red", xlab="1st. Digit", ylab = "Relative frequency",
pch=19, xaxp=c(1,9,8), ylim=c(0,0.35), main="Figure 4:
Relative Frequency of 1st. Digits of Coefficients,
Standard Errors, & t-Statistics",cex.main=0.9,cex.lab=0.9,cex.axis=0.8)
points(d1,bfds$bfd$data.dist, col="blue", pch=8)
points(d1,bfdt$bfd$data.dist, col="black", pch=17)
lines(d1, bfdc$bfd$benford.dist)
legend(6.5,0.35, legend=c("Benford's law", "coefficients","standard errors", "t-statistics"),
       col=c("black", "red", "blue","black"), lty=c(1,NA,NA,NA),pch=c(NA,19,8,17), cex=0.8)

#####################
#
#SECOND DIGITS
#

coef_2_ch <- str_sub(as.character(abs(coef_total)), 2, 2)    # Need to take abs. value to eliminate the "-" CHARACTERS
se_2_ch<- str_sub(as.character(se_total), 2, 2)
tstat_2_ch<- str_sub(as.character(abs(tstat_total)), 2, 2)
# Convert back to numeric
coef_2 <- as.numeric(coef_2_ch)
se_2<- as.numeric(se_2_ch)
tstat_2<- as.numeric(tstat_2_ch)
# Remove the "NAs"
coef_2 <- coef_2[!is.na(coef_2)]
se_2<- se_2[!is.na(se_2)]
tstat_2<- tstat_2[!is.na(tstat_2)]

coef_share<- c()
se_share<- c()
tstat_share<- c()

for (ii in 1:10) {                                  
i1<- ii-1                      # run from zero to 9
coef_share[ii] <- sum(coef_2 == i1)/length(coef_2)
se_share[ii] <- sum(se_2 == i1)/length(se_2)
tstat_share[ii] <- sum(tstat_2 == i1)/length(tstat_2)
}

plot(d2,coef_share, col="red", xlab="2nd. Digit", ylab = "Relative frequency",
pch=19, xaxp=c(0,9,9), ylim=c(0.08,0.14), main="Figure 5:
Relative Frequency of 2nd. Digits of Coefficients,
Standard Errors, & t-Statistics",cex.main=0.9,cex.lab=0.9,cex.axis=0.8)
lines(d2, benford_2)
points(d2,se_share, col="blue", pch=8)
points(d2,tstat_share, col="black", pch=17)
legend(6.5,0.13, legend=c("Benford's law", "coefficients","standard errors", "t-statistics"),
       col=c("black", "red", "blue","black"), lty=c(1,NA,NA,NA),pch=c(NA,19,8,17), cex=0.8)
#######################################################
