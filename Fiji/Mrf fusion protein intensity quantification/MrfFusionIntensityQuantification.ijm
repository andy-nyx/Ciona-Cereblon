//This script is a modification of the post of Thomas Boudier (https://forum.image.sc/t/3d-roi-manager-measurements-multi-channel/20164)

// use plugins/macros/record to record the necessary opening and channel splitting
// 3 images are opened in this example : labels.tif with segmented object, C1-signal.tif and C2-signaltif with signal to quantify
// select measurments in Manager3D Options, select Plugins/Record to record
dir1 = getDirectory("Choose Source Directory for image to be quantified"); //source image
index = lastIndexOf(dir1,"\\");
print(index);
print(dir1);
tempdir = substring(dir1,0,index);
dir2 = tempdir+"_CheckBinary\\"; //source binary
print(dir2);
dir3 = tempdir+"_quantification\\";
File.makeDirectory(dir3);
run("3D Manager Options", "volume integrated_density mean_grey_value std_dev_grey_value mode_grey_value minimum_grey_value maximum_grey_value distance_between_centers=10 distance_max_contact=1.80 drawing=Contour");// run the manager 3D and add image
run("3D Manager");
j = getFileList(dir1);
for (i=0; i<j.length; i++)
{
	open(dir1 + j[i]);
	title1=getTitle();
	print(title1);
	rename("signal.tif");
	index1 = lastIndexOf(title1, ".");
	title2 =  substring(title1,0,index1);
	open(dir2 + "NucleiBinary_"+title1);
	selectWindow("NucleiBinary_"+title1);
	rename("labels.tif");
// select the iamge with the labelled objects
	selectWindow("labels.tif");
	Ext.Manager3D_Segment(128, 255);
	Ext.Manager3D_AddImage();
	selectWindow("signal.tif");
	run("Split Channels");
// Quantifications, save measurements and close window
	selectWindow("C1-signal.tif"); // nG intensity
	Ext.Manager3D_Quantif();
	Ext.Manager3D_SaveResult("Q",dir3+"C1-nG-Quantif3D"+title2+".txt");
	Ext.Manager3D_CloseResult("Q"); 
	selectWindow("C2-signal.tif"); // mCh intensity
	Ext.Manager3D_Quantif();
	Ext.Manager3D_SaveResult("Q",dir3+"C2-mCh-Quantif3D"+title2+".txt");
	Ext.Manager3D_CloseResult("Q");
	Ext.Manager3D_Erase();
	run("Clear Results");
	run("Close All");
}