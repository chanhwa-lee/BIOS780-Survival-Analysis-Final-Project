data data;
	infile "/home/u49431079/final_data.txt" firstobs=2;
	input treat$ extloc$ regimp$ noposnd age time_dea status_dea log_noposnd;
run;



proc phreg data=data;
	class treat(ref="1") extloc(ref="1") regimp(ref="0");
	model time_dea*status_dea(0) = treat extloc regimp age extloc*age log_noposnd;
	assess var=(age log_noposnd) ph / resample;
	output out=res RESSCO=a b c d e f g h i j k l m n o p q r s t u v;
run;



proc sort data=res;
	by time_dea;
run;


title "Score residual";
proc sgplot data=res;
  series x=time_dea y=a;
run;
title;

title "Score residual";
proc sgplot data=res;
  series x=time_dea y=b;
run;
title;

title "Score residual";
proc sgplot data=res;
  series x=time_dea y=c;
run;
title;

title "Score residual";
proc sgplot data=res;
  series x=time_dea y=d;
run;
title;

title "Score residual";
proc sgplot data=res;
  series x=time_dea y=e;
run;
title;

title "Score residual";
proc sgplot data=res;
  series x=time_dea y=f;
run;
title;

title "Score residual";
proc sgplot data=res;
  series x=time_dea y=g;
run;
title;

title "Score residual";
proc sgplot data=res;
  series x=time_dea y=h;
run;
title;

title "Score residual";
proc sgplot data=res;
  series x=time_dea y=i;
run;
title;

title "Score residual";
proc sgplot data=res;
  series x=time_dea y=j;
run;
title;

title "Score residual";
proc sgplot data=res;
  series x=time_dea y=k;
run;
title;