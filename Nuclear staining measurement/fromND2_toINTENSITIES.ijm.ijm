// Batch convert Nikon nd2 files to tiff files under one selected dir,
// sub directories are ignored.

input = getDirectory("source directory");
output = getDirectory("output directory");
list = getFileList(input);

  for (i=0; i<list.length; i++) {
	image = list[i];
    path = input + image;
    print(path);
    
    run("Bio-Formats Macro Extensions");
    Ext.setId(path);
    Ext.getCurrentFile(file);
    //get how many series are in the nd2 file
    Ext.getSeriesCount(seriesCount);
	print("Series Count: "+seriesCount);

    for (s=1; s<=seriesCount; s++)                                                  
    {
        run("Bio-Formats Importer", "open=&path autoscale color_mode=Default view=Hyperstack stack_order=XYCZT series_"+s);
        t=getTitle();
        close();
        
        run("Bio-Formats Importer", "open=&path autoscale color_mode=Default split_channels view=Hyperstack stack_order=XYCZT series_"+s);  
        
		print(t);
        selectWindow(t + " - C=0");
        //setOption("ScaleConversions", true);
		//run("8-bit");
        run("Z Project...", "projection=[Max Intensity]");
		rename("histone");
		selectWindow(t + " - C=1");
		//setOption("ScaleConversions", true);
		//run("8-bit");
		run("Duplicate...", "title=nuclei");

		setAutoThreshold("Otsu dark");
		//run("Threshold...");
		setOption("BlackBackground", true);
		run("Convert to Mask");
		run("Watershed");

		run("Set Measurements...", "area mean min median limit display redirect=histone decimal=3");
		run("Analyze Particles...", "size=30-Infinity display exclude add");

        run("Close All");
    }
  }

saveAs("Results", output + "/Results.csv");