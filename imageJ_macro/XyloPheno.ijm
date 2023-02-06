/////////////////////
// Xylem phenology //
/////////////////////

////
// Developed by Estelle Noyer
// Version 7 _ October 2020
// Czechglobe institute, BRNO, Czechglobe
// noyer.estelle@gmail.com
////


// Image preparation // ------------------------------- //

// Selection
scale = getBoolean("Did you set the scale?")
if(scale == false){
	run("Set Scale...");	
}

setTool("polygon");
title = "Pause";
msg = "Delineate the tree-ring or the interested zone.\nTo select the entire picture: ctl + a";
waitForUser(title, msg);

// Duplicate
run("Duplicate...", " ");
waitForUser("Save (ctrl+s)");
name_image = getInfo("image.filename");
pathway = getInfo("image.directory");
roiManager("Add");
roiManager("Select", 0);
roiManager("Rename", "Zone");
roiManager("Save", pathway + name_image + ".zip");

// Export area of the zone
run("Set Measurements...", "area perimeter bounding feret's invert redirect=None decimal=3");
roiManager("Measure");
saveAs("Results", pathway + name_image + "_zone.csv");

// preparation ending
run("Clear Outside");
run("Select None");

roiManager("reset");
run("Clear Results");

// Cells counting // ------------------------------- //
//first round
Ccells = getBoolean("Do you want to count cells ?");

if(Ccells == true){
	// Cambial Row //
		// point tool
	do{
		setTool("point");
		run("Point Tool...", "type=Dot color=Yellow size=Small add label");
		roiManager("Show All");
		title = "Pause";
		msg = "Select cambial rows\nor deselect cells with Alt+clicking on the dot";
		waitForUser(title, msg);
		check = getBoolean("Are you satisfied?");
	}while(check == false);
	
		// name ROI
	labelRM = "CR";
	nb_initial = roiManager("count");
	for(i=0; i<nb_initial; i++){
		roiManager("Select", i);
		roiManager("Rename", labelRM + "." + i+1);
		run("Set Measurements...", "centroid invert redirect=None decimal=3");
		run("Measure");
	}
	roiManager("Save", pathway + name_image + "_cellcount.zip");

		// name results  + count + percentage
	labelC = "Cambial Row";
	row_ini = nResults;
	for(i=0; i<row_ini; i++){
		setResult("Cells", i, labelC);
	}
	saveAs("Results", pathway + name_image + "_cellcount.csv");

	// Xylem elements //
	labelRM = newArray("APF.E", "V.E", "APF.SCW", "V.SCW", "F.M", "AP.M", "V.M");
	labelC = newArray("01_Enlarging Axial Parenchyma_Fibers", "02_Enlarging Vessels", "03_SCW Axial Parenchyma_Fibers", "04_SCW Vessels", "05_Mature Fibers", "06_Mature Axial Parenchyma", "07_Mature Vessels");

	for(j=0; j<lengthOf(labelRM); j++){
		
		cell = getBoolean("Do you want to count " + labelC[j] + "?");
		
		if(cell == true){
			// point tool
			do{
				setTool("point");
				run("Point Tool...", "type=Dot color=Green size=Small add label");
				roiManager("Show All");
				title = "Pause";
				msg = "Select " + labelC[j] + "\nor deselect cells with Alt+clicking on the dot";
				waitForUser(title, msg);
				check = getBoolean("Are you satisfied?");
			}while(check == false);
				
			// name ROI
			tot = roiManager("count");
			nb = tot - nb_initial;
			for(i=nb_initial; i<tot; i++){
				roiManager("Select", i)
				roiManager("Rename", labelRM[j] + "." + i+1);
				run("Set Measurements...", "centroid invert redirect=None decimal=3");
				run("Measure");
			}
			roiManager("Save", pathway + name_image + "_cellcount.zip");
			nb_initial = roiManager("count");
			
			// name results + count
			row = nResults;
			for(i=row_ini; i<row; i++){
				setResult("Cells", i, labelC[j]);
			}
			saveAs("Results", pathway + name_image + "_cellcount.csv");
			row_ini = nResults;
		}
	}

}else{}


//second round
VerifC = getBoolean("Do you want to make a second round ?");

if(VerifC == true){
	labelRM = newArray("APF.E", "V.E", "APF.SCW", "V.SCW", "F.M", "AP.M", "V.M");
	labelC = newArray("01_Enlarging Axial Parenchyma_Fibers", "02_Enlarging Vessels", "03_SCW Axial Parenchyma_Fibers", "04_SCW Vessels", "05_Mature Fibers", "06_Mature Axial Parenchyma", "07_Mature Vessels");
	
	for(j=0; j<lengthOf(labelRM); j++){
		
		cell = getBoolean("Do you want to count " + labelC[j] + "?");
		
		if(cell == true){
			// point tool
			do{
				setTool("point");
				run("Point Tool...", "type=Dot color=Green size=Small add label");
				roiManager("Show All");
				title = "Pause";
				msg = "Select " + labelC[j] + "\nor deselect cells with Alt+clicking on the dot";
				waitForUser(title, msg);
				check = getBoolean("Are you satisfied?");
			}while(check == false);
				
			// name ROI
			tot = roiManager("count");
			nb = tot - nb_initial;
			for(i=nb_initial; i<tot; i++){
				roiManager("Select", i)
				roiManager("Rename", labelRM[j] + "." + i+1);
				run("Set Measurements...", "centroid invert redirect=None decimal=3");
				run("Measure");
			}
			roiManager("Save", pathway + name_image + "_cellcount.zip");
			nb_initial = roiManager("count");
			
			// name results + count
			row = nResults;
			for(i=row_ini; i<row; i++){
				setResult("Cells", i, labelC[j]);
			}
			saveAs("Results", pathway + name_image + "_cellcount.csv");
			row_ini = nResults;
		}
	}
	
	roiManager("reset");
	run("Clear Results");
	
}else{
	roiManager("reset");
	run("Clear Results");}	


// Vessels // ------------------------------- //
vessels = getBoolean("Do you want to measure vessel size?");

if(vessels == true){
	// image preparation
	selectImage(name_image);
	run("Duplicate...", " ");
	run("16-bit");
	Vessel_image = name_image + "_Vessel";
	
	saveAs("tiff", pathway + Vessel_image + ".tif");

	run("Threshold...");
	title = "Pause";
	msg = "Adjust the threshold then tap on OK";
	waitForUser(title, msg);
	
	roiManager("Open", pathway + name_image + ".zip");
	roiManager("Save", pathway + name_image + "_Vessel.zip");

	// lumen area
	do {
		roiManager("reset");
		roiManager("Open", pathway + name_image + "_Vessel.zip");
		roiManager("Show All");
		run("Select None");
		setTool("wand");
		title = "Pause";
		msg = "Select all vessels \nAdd (selecting an area then click on [t]) or \nErase (selecting an area then click on [Delete])";
		waitForUser(title, msg);
		roiManager("Save", pathway + name_image + "_Vessel.zip");
		vessel_ok = getBoolean("satisfait?");	
	} while (vessel_ok == false);
	roiManager("Save", pathway + name_image + "_Vessel.zip");
	roiManager("Open", pathway + name_image + ".zip");
	run("Set Measurements...", "area centroid feret's invert redirect=None decimal=3");
	roiManager("Measure");
	saveAs("Results", pathway + name_image + "_Vessel.csv");

	nbV = roiManager("count") - 1;

	roiManager("reset");
	run("Clear Results");
	

//////////////////////////////////////////////////////////////////////////// END /////////