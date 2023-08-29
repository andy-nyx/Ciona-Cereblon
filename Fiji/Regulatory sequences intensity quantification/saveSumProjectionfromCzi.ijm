//All the images that you want to quantify are open in fiji

dir1 = getDirectory("Choose Source Directory");
index = lastIndexOf(dir1,"\\");
tempdir = substring(dir1,0,index);
dir2 = tempdir+"/RegulatorySequences_sum/"
File.makeDirectory(dir2);
j=nImages;
print(j);
for (i=0; i<j; i++)
{
title1=getTitle();
	index2 = lengthOf(title1);
	print(index2);
	index3 = lastIndexOf(title1, " - ");
	print(index3);
	title2 = substring(title1,index3+3,index2);
run("Z Project...", "projection=[Sum Slices]");
saveAs("tiff", dir2+"Sum_"+title2+".tif");
close();
close(title1);
}
