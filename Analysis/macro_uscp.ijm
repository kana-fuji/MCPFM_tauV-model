//
//run: imagej -b macro_usc.ijm 00 
//since:20200324 update:20200717
//Kana Fuji
//
//main------------------------------------------------------start
param = getArgument;
if (param=="") exit ("No argument!");

inputdir = "../IN/"+param;
datadir = "../DATA/"+param;
outputdirpng = "PNG";
outputdiravi = "AVI";
outputdirgif = "GIF";

sign=uscCheck(datadir);
if (sign==0) exit ("No data!");
//print(sign);

tmax=loadTime(datadir);
tmin=0;

//tout=loadParam(inputdir);
di=loadParam(datadir);
//di=tout;
//print(di);

//read time series----------------------start

for (i=tmin; i<=tmax; i=i+di) {
istr = "" + i;
run("Text Image... ", "open=[" + datadir + "/u_" + istr + ".dat]");
}
run("Images to Stack", "name=cell title=[] use");
setMinAndMax(0.2000, 1.0000);

if(sign==2 || sign==3 || sign==5){
for (i=tmin; i<=tmax; i=i+di) {
istr = "" + i;
run("Text Image... ", "open=[" + datadir + "/s_" + istr + ".dat]");
}
run("Images to Stack", "name=lumen title=[] use");
setMinAndMax(0.5000, 1.0000);
}

if(sign==3 || sign==4 || sign==5){
for (i=tmin; i<=tmax; i=i+di) {
istr = "" + i;
run("Text Image... ", "open=[" + datadir + "/c_" + istr + ".dat]");
}
run("Images to Stack", "name=ecm title=[] use");
setMinAndMax(0.5000, 1.0000);
}

if(sign==5){
for (i=tmin; i<=tmax; i=i+di) {
istr = "" + i;
run("Text Image... ", "open=[" + datadir + "/p_" + istr + ".dat]");
}
run("Images to Stack", "name=anti title=[] use");
setMinAndMax(0.5000, 1.0000);
}

//read time series--------------------end


//merge channels----------------------start

if(sign==5){
run("Merge Channels...", "c1=[cell]  c2=[anti]   c5=[lumen]   c7=[ecm]   create");
}

if(sign==3){
run("Merge Channels...", "c1=[cell]   c5=[lumen]   c7=[ecm]   create");
}
else if(sign==2){
run("Merge Channels...", "c1=[cell]   c5=[lumen]   create");
}
else if(sign==4){
run("Merge Channels...", "c1=[cell]   c7=[ecm]   create");
}
else if(sign==1){
run("Red");
}

//merge channels----------------------end

//timestamp---------------------------start
Stack.getDimensions(w, h, channels, slices, frames);
//Stack.getPosition(ch, sl, fr);
//dt=tout[0]*tout[1];
//print(dt);

if(sign==3||sign==5){setForegroundColor(0, 255, 0);}
//if(sign==1){setForegroundColor(0, 0, 0);}
//run("Label...", "format=0000 starting=0.00 interval=s x=30 y=75 font=30 text=[] range=1-slices");
para="format=0000 starting=0.0 interval="+di+" x=30 y=75 font=30 text=[] range=1-"+slices;
run("Label...", para);



//for (i=tmin; i<=tmax; i=i+di) {
//text="t="+i*tout[0];
//makeText(text,10,20);
//setFont("Century", 60, " antialiased");
//run("Add Selection...", "stroke=black fill=yellow new");
//run("Select None");
//j++;
//}

//timestamp---------------------------end

//output image------------------------start

saveAs("AVI", outputdiravi+"/"+param+".avi");
//saveAs("GIF", outputdirgif+"/"+param+".gif");
setSlice(nSlices);
saveAs("PNG", outputdirpng+"/"+param+"_ij.png");
//saveAs("Tiff", outputdirectory + "/" + prefix + filelist[i]);

close();

//main------------------------------------------------------------------end


//function--------------------------------------------------------------start

function loadTime(data){
time=0;
//print(data);
if(File.exists(data+"/timestamp.dat")){
str = File.openAsString(data+"/timestamp.dat");
//print(str);
parLines=split(str,"\n");
time=parseInt(parLines[0]);
}
 return time;
}

function loadParam(data){
tout=0;
if(File.exists(data+"/param.txt")){
str = File.openAsString(data+"/param.txt");
parLines=split(str,"\n");
tout=parseInt(parLines[0]);
}
 return tout;
}

function uscCheck(dir){
v=0;
if(File.exists(dir+"/c_0.dat")){
	if(File.exists(dir+"/s_0.dat")){
	v=3;
		if(File.exists(dir+"/p_0.dat")){
			v=5;
		}
	}
	else{
	v=4;
	}
}
else if(File.exists(dir+"/s_0.dat")){
v=2;
}
else{
v=1;
}
 return v;
}
