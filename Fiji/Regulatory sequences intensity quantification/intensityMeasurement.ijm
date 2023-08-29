dir1 = getDirectory("Choose Source Directory for Max projection"); //max projection
index = lastIndexOf(dir1,"\\");
print(index);
tempdir = substring(dir1,0,index);
dir2 = tempdir+"_mask/"
dir3 = getDirectory("Choose Source Directory for Sum projection"); //sum projection
print(dir2);
print(dir1);
print(dir3);
File.makeDirectory(dir2);
j = getFileList(dir1);
l = getFileList(dir3);
for (i=0; i<j.length; i++)
	{open(dir1 + j[i]);
	title1=getTitle();

	index1 = lastIndexOf(title1, ".");
	title2 =  substring(title1,0,index1);
	print(title2);
	selectWindow(title1);
	title = "Select area";
	msg = "Select area";
	waitForUser(title, msg);
	run("Clear Outside");
	run("Split Channels");
	close("C2-"+title2+".tif");
	selectWindow("C1-"+title2+".tif");
	selectWindow("C1-"+title2+".tif");
	run("Duplicate...", " ");
	selectWindow("C1-"+title2+"-1.tif");
	run("Median...", "radius=30");
	imageCalculator("Subtract create", "C1-"+title2+".tif","C1-"+title2+"-1.tif");
	selectWindow("Result of C1-"+title2+".tif");
	run("Smooth");
	run("Smooth");
	run("Smooth");
	run("Smooth");
	setMinAndMax(0, 65535); 
	setOption("ScaleConversions", true);
	run("8-bit");
	run("Auto Local Threshold", "method=Bernsen radius=15 parameter_1=0 parameter_2=0 white");
	run("Create Selection");
	run("Create Mask");
	saveAs("tiff", dir2+"mask_"+title2+".tif");
	open(dir3 + l[i]);
	title3 = getTitle();
	run("Split Channels");
	close("C2-"+title3);
	close("C3-"+title3);
	selectWindow("C1-"+title3);
	run("Restore Selection");
	run("Measure");
	run("Close All");
	}
saveAs("Results", dir3+"Results.csv");
run("Clear Results");
k = getFileList(dir2);
for (i=0; i<k.length; i++)
{open(dir2 + k[i]);
	title1=getTitle();
	run("Watershed");
	run("Analyze Particles...", "size=10-75 circularity=0.1-1.00 display summarize");
	run("Close All");
}
selectWindow("Summary");
saveAs("Results", dir2+"nucleusSummary.csv");
selectWindow("Results");
saveAs("Results", dir2+"nucleusResults.csv");