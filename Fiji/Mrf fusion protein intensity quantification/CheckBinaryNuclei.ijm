dir1 = getDirectory("Image Folder"); //Images used to make the binary
dir2 = getDirectory("Binary Folder"); //Binary to check
index = lastIndexOf(dir1,"\\");
tempdir = substring(dir1,0,index);
dir3 = tempdir+"_CheckBinary/";
File.makeDirectory(dir3);
j = getFileList(dir1);
for (i=0; i<j.length; i++)
{
	open(dir1 + j[i]);
	title1=getTitle();
	index1 = lastIndexOf(title1, ".");
	title2 =  substring(title1,0,index1);
	print(title2);
	open(dir2 + "NucleiBinary_"+title1);
	run("16-bit");
	run("Invert LUT");
	rename("NucleiBinary");
	selectWindow(title1);
	run("Split Channels");
	selectWindow("C2-"+title1);
	rename("Nuclei");
	run("Merge Channels...", "c1=NucleiBinary c2=Nuclei create keep");
	Stack.setChannel(1);
	setMinAndMax(0, 500);
	Stack.setChannel(2);
	setMinAndMax(0, 25000);
	waitForUser("Check Binary");
	run("Split Channels");
	selectWindow("C1-Composite");
	setThreshold(128, 65535, "raw");
	run("Make Binary", "method=Default background=Default create");
	saveAs("tiff", dir3+"NucleiBinary_"+title2+".tif");
	run("Close All");
}	