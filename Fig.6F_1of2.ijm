// This script cuts out single cells in 3D and easure nascent RAN-FISH spot intensity   /23-09-2020/
// input 1: where directory with FOVs.tif are deposited in single folders (use seperate script for that!)
// input 2: standard singlecellmask folder with subdirs
// output: single cell images of smFISH, cropped to 500x500px (Saves space!)





//1. DEFINE FUNCTION PRODUCING 3D SINGLE CELL STACKS //
///////////////////////////////////////////////////////

end_TRITC=".tif";
end_Masks="_SingleCellMask.tif";

function multiplier(inputTRITC, inputMasks) {
setBatchMode(true);
listTRITC = getFileList(inputTRITC);
listMaskFiles = getFileList(inputMasks);
for (i = 0; i < listTRITC.length; i++){     
				if (endsWith(listTRITC[i], ".tif")){
					core = substring(listTRITC[i],0,lengthOf(listTRITC[i])-4);    // name ending: ".tif"
					//print(listTRITC.length);
					open(inputTRITC + core + end_TRITC);
					
							run("Make Substack...", "  slices=91-120");      // here change for INT-FISH Z-slices!!!!!!!
							run("Subtract Background...", "rolling=4 stack");
							run("Properties...", "channels=1 slices=30 frames=1 unit=px pixel_width=1.0000 pixel_height=1.0000 voxel_depth=1.0000"); // new 20-09-2020
							imgName_TRITC=getTitle();
				}
}

for (k = 0; k < listMaskFiles.length; k++){         // k - single cell masks
				if (endsWith(listMaskFiles[k], "_SingleCellMask.tif")){
						core = substring(listMaskFiles[k],0,lengthOf(listMaskFiles[k])-19); 
   						open(inputMasks + core + end_Masks);
						run("Select All");
						run("Invert");
						run("Multiply...", "value=0.0039");
						imgName_Mask=getTitle();  					

   					    imageCalculator("Multiply create stack", imgName_TRITC ,imgName_Mask);
   					imgName_Final=getTitle();		

                     //Cropper+Aligner below:
							run("Clear Results");
      					    title=getTitle;selectWindow(title);nz=nSlices;
      					    run("Set Measurements...", "center redirect=None decimal=4");
      					    run("Properties...", "channels=1 slices=1 frames=30 unit=px pixel_width=1 pixel_height=1 voxel_depth=1"); // note frames here changed 1->30
      					    run("Measure");
								x=getResult("XM",i-1); //translation vectors!
              				    y=getResult("YM",i-1);
               				    if(i==1){x0=x;y0=y;}
              				    sx=x0-x;sy=y0-y;
               				    run("Translate...", "x="+sx+" y="+sy+" slice");

							    run("Measure");  //measures center of gravity in aligned image
      								 x1=getResult("XM",i-1); //translation vectors!
      							     y1=getResult("YM",i-1);
         							 print(x1);
         							 print(y1);
        							 setTool("rectangle");
									 run("Specify...", "width=220 height=220 x=x1 y=y1 slice=10 constrain centered scaled"); 
									 run("Crop"); //crops the image to a single cell

               				    	 //selectWindow(title);
               				   		 saveAs("Tiff", inputMasks+core+"_smFISH_INT");	

		
				}
}
setBatchMode(false);
close();
} //end of function----------------------------------------------------------------------



//2. PERFORM CUTTING OUT SINGLE 3D CELLS AND SAVING THEM //
///////////////////////////////////////////////////////////


// 0. ADD A STEP WHERE TIFF IMAGES ARE DEPOSITED INTO FOLDERS!!!!!!!!!!!!!! //


mainDirTRITC = getDirectory("Choose a main directory "); //main directory containing single file folders
mainTRITCList = getFileList(mainDirTRITC);
print(mainTRITCList.length); //to make sure it recognized folders

mainDirMasks = getDirectory("Choose a main directory "); //main directory containing single file folders
mainMasksList = getFileList(mainDirMasks);
print(mainMasksList.length); //to make sure it recognized folders

setBatchMode(true);
for (l=0; l<mainTRITCList.length; l++) {  // for loop to parse through names in main folder
     if(endsWith(mainTRITCList[l], "/")){   // if the name is a subfolder...

          subDir = mainDirTRITC + mainTRITCList[l]; //directory of l-folder, - one of the inputs of multiplier()
          subDirMasks = mainDirMasks + mainMasksList[l];
          multiplier(subDir, subDirMasks);

     }
}
setBatchMode(false);
wait(500);



//3. PERFORM 3D-OC FOR INT-SPOT INTENSITY //
////////////////////////////////////////////
setBatchMode(true);
for (b=0; b<mainMasksList.length; b++) { 
singleDir=mainDirMasks + mainMasksList[b];  //define a folder directory
listInDir = getFileList(singleDir);
print(listInDir.length);
for (d=0; d<listInDir.length; d++) { 
if (endsWith(listInDir[d], "INT.tif")){
	print(listInDir[d]);
	open(singleDir+listInDir[d]); 
	cell=getTitle;
	run("Properties...", "channels=1 slices=30 frames=1 unit=px pixel_width=1.0000 pixel_height=1.0000 voxel_depth=1.0000");
	//getStatistics(area, mean, min, max, std, histogram);
	Stack.getStatistics(count, mean, min, max, std);      // get statistics from the ENTIRE stack, not only the first frame!
	if (max < 105) threshold = max;                       // if threshold is the max of an image it wont pick any spots anyways!
  	else threshold = 105;
  	print(threshold);
  	    // below tested settings, min 20-150, treshold 105 - works robust, max 150px^3 size might not be enough for larger bursts, then use 250!
  		run("3D Objects Counter", "threshold=" + threshold +" slice=13 min.=20 max.=190 exclude_objects_on_edges statistics"); //dont exclude objects on edges!	min 30 better!
		g=b+1; //numerator of FOVs
		//saveAs("txt",input+i+"_frame_"+x+"min"); close("statistics"); //close(cell); // TXT
		saveAs("Results",singleDir+listInDir[d]+g+"_fov_" + d + "_cell_3DNascentSpots.csv"); close("statistics"); //close(cell); // CSV
		
        // Do OC for 2D MaxProj image as it is the most precise in identification of the real nascent spots:
		selectWindow(cell);
		run("Z Project...", "projection=[Max Intensity]");
		run("Subtract Background...", "rolling=20");
		cell2D=getTitle;
		getStatistics(area, mean, min, max, std, histogram);
		if (max < 100) threshold2 = max; 
  		else threshold2 = 100;
		run("3D Objects Counter", "threshold=" + threshold2 +" slice=1 min.=11 max.=4000 exclude_objects_on_edges statistics"); 
		saveAs("Results",singleDir+listInDir[d]+g+"_fov_" + d + "_cell_2DNascentSpots.csv"); close("statistics"); close(cell); close(cell2D);

		//close immediately 3D-OC result (slow):
        	listWin = getList("window.titles"); //setBatchMode(true);
        	for (j=0; j<listWin.length; j++){
        	winame = listWin[j];
        	selectWindow(winame);
        	run("Close");
        	} //setBatchMode(false);

}
}
}
setBatchMode(false);

		   //close all 3D-OC results left:
        	listWin = getList("window.titles"); setBatchMode(true);
        	for (j=0; j<listWin.length; j++){
        	winame = listWin[j];
        	selectWindow(winame);
        	run("Close");
        	} setBatchMode(false);