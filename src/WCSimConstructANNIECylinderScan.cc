//	-*- mode:c++; tab-width:4;	-*-
#include "WCSimDetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"
#include "G4Torus.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4IntersectionSolid.hh"
#include "G4Polyhedra.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4PVReplica.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4AssemblyVolume.hh"
#include "G4SubtractionSolid.hh"
#include "globals.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpticalSurface.hh"
#include "G4UserLimits.hh"
#include "G4ReflectionFactory.hh"
#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "WCSimTuningParameters.hh" //jl145

#include <fstream>

/***********************************************************
 *
 * This file containts the functions which construct a 
 * cylindrical WC detector.	It used by both the SK and 
 * LBNE WC detector modes.	It is called in the Construct()
 * method in WCSimDetectorConstruction.cc.
 *
 * Sourcefile for the WCSimDetectorConstruction class
 *
 * This modified file (ConstructANNIECylinderScan) is designed 
 * explicitly for constructing the ANNIE detector, with 
 * the PMT positions extracted from laser scan data. 
 ***********************************************************/

G4LogicalVolume* WCSimDetectorConstruction::ConstructANNIECylinderScan()
{
	G4cout << "**** Building Cylindrical Detector ****" << G4endl;
	G4cout << "ConstructANNIECylinderScan" << G4endl;

	debugMode = false;
	
	//-----------------------------------------------------
	//---------------------Steel Barrel-------------------- 
	//-----------------------------------------------------
	
	totalAngle = 2.0*pi*rad;
	dPhi = totalAngle/ WCBarrelRingNPhi;
	
	WCIDRadius = WCIDDiameter/2.;
	WCLength = WCIDHeight;
	
	G4double WChalfheightexcess, WCradialexcess;
	// tank is 7GA (7-gauge steel) = (0.1793"normal sheet steel or) 0.1875" stainless steel
	// ^ gauge represents the thickness of 1 ounce of copper rolled out to an area of 1 square foot
	// = 4.7625mm thick 'A36' steel (density 7,800 kg/m3 )
	WChalfheightexcess=4.7625*mm; 
	WCradialexcess=4.7625*mm;
	
	G4Tubs* solidWC = new G4Tubs("WC",
								 0.0*m,
								 WCRadius+WCradialexcess, 
								 0.5*WCLength+WChalfheightexcess,
								 0.*deg,
								 360.*deg);
	
	G4LogicalVolume* logicWC;
	logicWC = new G4LogicalVolume(solidWC,
								  G4Material::GetMaterial("StainlessSteel"),
								  "WC",
								  0,0,0);
	
	G4VisAttributes* showColor = new G4VisAttributes(G4Colour(0.0,1.0,0.0));
	logicWC->SetVisAttributes(showColor);
	
	//logicWC->SetVisAttributes(G4VisAttributes::Invisible); //amb79
	
	//-----------------------------------------------------
	//---------------Water in the Barrel-------------------
	//-----------------------------------------------------
	
	// Decide if adding Gd
	water = (WCAddGd) ? "Doped Water" : "Water";
	
	G4Tubs* solidWCBarrel = new G4Tubs("WCBarrel",
										0.0*m,
										WCRadius,
										0.5*WCLength,
										0.*deg,
										360.*deg);
	
	logicWCBarrel = 
		new G4LogicalVolume(solidWCBarrel,
			G4Material::GetMaterial(water),
			"WCBarrel",
			0,0,0);
	
	physiWCBarrel = 
		new G4PVPlacement(0,
						  G4ThreeVector(0.,0.,0.),
						  logicWCBarrel,
						  "WCBarrel",
						  logicWC,
						  false,
						  0);
	
	// This volume needs to made invisible to view the blacksheet and PMTs with RayTracer
	if (Vis_Choice == "RayTracer"){
		logicWCBarrel->SetVisAttributes(G4VisAttributes::Invisible);
	} else {
		//{if(!debugMode)
		//logicWCBarrel->SetVisAttributes(G4VisAttributes::Invisible);} 
	}
	
	//-----------------------------------------------------------
	// The Blacksheet, a daughter of the cells containing PMTs,
	// and also some other volumes to make the edges light tight
	//-----------------------------------------------------------
	
	G4double annulusBlackSheetRmax[2] = {(WCIDRadius+WCBlackSheetThickness),
										  WCIDRadius+WCBlackSheetThickness};
	G4double annulusBlackSheetRmin[2] = {(WCIDRadius),
										  WCIDRadius};
	mainAnnulusHeight = WCIDHeight -2.*WCBarrelPMTOffset;
	G4double mainAnnulusZ[2] = {-mainAnnulusHeight/2., mainAnnulusHeight/2};
	
	G4Polyhedra* solidWCBarrelCellBlackSheet = new G4Polyhedra("WCBarrelCellBlackSheet",
																WCBarrelCellStartPhi, // phi start
																totalAngle,           //total phi
																WCBarrelRingNPhi,     //NPhi-gon
																2,
																mainAnnulusZ,
																annulusBlackSheetRmin,
																annulusBlackSheetRmax);
	
	logicWCBarrelCellBlackSheet =
		new G4LogicalVolume(solidWCBarrelCellBlackSheet,
							G4Material::GetMaterial("Blacksheet"),
							"WCBarrelCellBlackSheet",
							0,0,0);
	
	G4VPhysicalVolume* physiWCBarrelCellBlackSheet =
		new G4PVPlacement(0,
						  G4ThreeVector(0.,0.,InnerStructureCentreOffset),
						  logicWCBarrelCellBlackSheet,
						  "WCBarrelCellBlackSheet",
						  logicWCBarrel,
						  false,
						  0,true);
	G4cout<<"Constructed barrel cell blacksheet with radius "<<WCIDRadius<<" to "<<(WCIDRadius+WCBlackSheetThickness)<<G4endl;
	WCCylInfo[2]=WCIDRadius/10.;
	
	G4LogicalBorderSurface * WaterBSBarrelCellSurface 
		= new G4LogicalBorderSurface("WaterBSBarrelCellSurface",
									 physiWCBarrel,
									 physiWCBarrelCellBlackSheet,
									 OpWaterBSSurface);
	
	// Change made here to have the if statement contain the !debugmode to be consistent
	// This code gives the Blacksheet its color. 
	
	if(Vis_Choice == "RayTracer"){
		G4VisAttributes* WCBarrelBlackSheetCellVisAtt 
			= new G4VisAttributes(G4Colour(0.2,0.9,0.2)); // green color
		WCBarrelBlackSheetCellVisAtt->SetForceSolid(true); // force the object to be visualized with a surface
		WCBarrelBlackSheetCellVisAtt->SetForceAuxEdgeVisible(true); // force auxiliary edges to be shown
		if(!debugMode)
			logicWCBarrelCellBlackSheet->SetVisAttributes(WCBarrelBlackSheetCellVisAtt);
		else
			logicWCBarrelCellBlackSheet->SetVisAttributes(G4VisAttributes::Invisible);
	} else {
		G4VisAttributes* WCBarrelBlackSheetCellVisAtt 
			= new G4VisAttributes(G4Colour(0.2,0.9,0.2));
		if(!debugMode)
			logicWCBarrelCellBlackSheet->SetVisAttributes(G4VisAttributes::Invisible);
		else
			logicWCBarrelCellBlackSheet->SetVisAttributes(WCBarrelBlackSheetCellVisAtt);
	}
	
	//-----------------------------------------------------
	//--------------------- The PMT -----------------------
	//-----------------------------------------------------
	
	G4LogicalVolume* logicWCPMT;
	std::vector<G4LogicalVolume*> logicWCPMTs;
	for(auto atankcollection : WCTankCollectionNames){
		G4String thepmtname = WCPMTNameMap.at(atankcollection);
		logicWCPMT= ConstructPMT(thepmtname, atankcollection, "tank");
		logicWCPMTs.push_back(logicWCPMT);
	}
	G4LogicalVolume* logicWCLAPPD=0;
	logicWCLAPPD = ConstructLAPPD(WCLAPPDName, WCIDCollectionName2);
	
	G4VisAttributes* WClogic 
			= new G4VisAttributes(G4Colour(0.4,0.0,0.8));
	WClogic->SetForceSolid(true);
	WClogic->SetForceAuxEdgeVisible(true);
	
	//logicWCPMT->SetVisAttributes(WClogic);
	//TODO: Re-enable Invisible attribute for PMTs after testing
	for(auto alogicWCPMT : logicWCPMTs) alogicWCPMT->SetVisAttributes(G4VisAttributes::Invisible);
	if(logicWCLAPPD) logicWCLAPPD->SetVisAttributes(G4VisAttributes::Invisible);
	
	//////--------------------------------------------------------------
	//////----------- Barrel + Top + Bottom PMT placement --------------
	//////--------------------------------------------------------------

	// Rotate in a way that the y axis points upwards (into the sky)
	G4RotationMatrix* WCPMTRotation = new G4RotationMatrix;
	WCPMTRotation->rotateY(90.*deg);
	
	// Create rotation matrices for the orientations of the PMTs
	std::vector<G4RotationMatrix*> pmt_rotation_matrices;
	// Bottom PMTs have panel number 0
	G4RotationMatrix *WCBottomCapRotation = new G4RotationMatrix();
	pmt_rotation_matrices.push_back(WCBottomCapRotation);	
	// Barrel PMTs have panel numbers 1-8
	for(int facei=0; facei<WCBarrelRingNPhi; facei++){
		G4RotationMatrix* WCPMTRotationNext = new G4RotationMatrix(*WCPMTRotation);
		WCPMTRotationNext->rotateX((dPhi*facei)-67.5*deg+180*deg);
		pmt_rotation_matrices.push_back(WCPMTRotationNext);
	}
	// Top PMTs have panel number 9
	G4RotationMatrix *WCTopCapRotation = new G4RotationMatrix();
	WCTopCapRotation->rotateY(180.*deg);
	pmt_rotation_matrices.push_back(WCTopCapRotation);

	G4cout <<"Size of pmt_rotation_matrices: "<<pmt_rotation_matrices.size()<<G4endl;

	// Read in file with PMT positions, create PMTs
	std::ifstream pmt_position_file("PMTPositions_Scan.txt");
	std::string next_pmt;
	double pmt_x, pmt_y, pmt_z, pmt_dirx, pmt_diry, pmt_dirz;
	double pmt_x_shift, pmt_y_shift, pmt_z_shift;
	int panel_nr, pmt_type;
	int PMTID;
	while (!pmt_position_file.eof()){
		pmt_position_file >> PMTID >> panel_nr >> pmt_x >> pmt_y >> pmt_z >> pmt_dirx >> pmt_diry >> pmt_dirz >> pmt_type;
		if (pmt_position_file.eof()) break;
		G4cout << "Read in PMT "<<PMTID<<", panel nr: "<<panel_nr<<", Position ("<<pmt_x<<","<<pmt_y<<","<<pmt_z<<"), PMT type: "<<pmt_type<<G4endl;
		G4LogicalVolume *logicWCPMT = logicWCPMTs.at(pmt_type);
		G4RotationMatrix *pmt_rot = pmt_rotation_matrices.at(panel_nr);
		pmt_x_shift = pmt_x*cm;
		pmt_y_shift = (168.1-pmt_z)*cm;
		pmt_z_shift = ((pmt_y+14.45))*cm;
		//pmt_z_shift = ((pmt_y+14.45)-InnerStructureCentreOffset/10.)*cm;
		G4cout <<"Edited PMT position ("<<pmt_x_shift<<","<<pmt_y_shift<<","<<pmt_z_shift<<")"<<G4endl;
		G4ThreeVector PMTPosition(pmt_x_shift,pmt_y_shift,pmt_z_shift);
		G4VPhysicalVolume *physicalWCPMT = new G4PVPlacement(pmt_rot,	//its rotation
															PMTPosition,		//its position
															logicWCPMT,			//its logical volume
															"WCPMT",			//its name
															logicWCBarrel,		//its mother volume
															false,				//no boolean operations
															PMTID,				//ID for this PMT (=channelkey in data)
															true);				//check overlaps*/
	}
	pmt_position_file.close();

	// Leave LAPPD placement in for now (?)
	
	//////-------------------------------------------------
	//////----------- Barrel LAPPD placement --------------
	//////-------------------------------------------------

	G4RotationMatrix* WCLAPPDRotation = new G4RotationMatrix;
	WCLAPPDRotation->rotateY(90.*deg);
	
	G4cout <<"Barrel LAPPD placement"<<G4endl;
	for(int facei=0; facei<WCBarrelRingNPhi; facei++){
		// rotate the LAPPD to point to tank centre. We need a new rotation matrix for each G4VPhysicalVolume
		G4RotationMatrix* WCLAPPDRotationNext = new G4RotationMatrix(*WCLAPPDRotation);
		WCLAPPDRotationNext->rotateX((dPhi*facei)-90.*deg);
		
		// position LAPPD on corner of the inner strucutre
		// Note PMTs use double facei=0.5, using int facei=0 accounts for the relative rotation
		G4double CellCentreX = WCIDRadius * sin(dPhi*facei);
		G4double CellCentreY = WCIDRadius * cos(dPhi*facei);
		
		double verticalSpacingLAPPD	= mainAnnulusHeight/(WCLAPPDperCellVertical+1);
		
		for(G4double j = 0; j < WCLAPPDperCellVertical; j++){	// num LAPPD cols in the central ring
		
		G4ThreeVector LAPPDPosition = G4ThreeVector(CellCentreX,
													CellCentreY,
													-mainAnnulusHeight/2.+(j+1.)*verticalSpacingLAPPD);
		
		G4VPhysicalVolume* physiWCBarrelLAPPD =
		new G4PVPlacement(WCLAPPDRotationNext,                      // its rotation
							LAPPDPosition,                          // its position
							logicWCLAPPD,                           // its logical volume
							"WCLAPPD",                              // its name
							logicWCBarrel,                          // its mother volume
							false,                                  // no boolean operations
							(int)(j+facei*WCLAPPDperCellVertical),  // a unique copy number
							false);                                 // don't check overlaps
		}
	}
	
	//--------------------------------------------------------------
	//---------------------------- CAPS ----------------------------
	//--------------------------------------------------------------

	// Leave the ConstructANNIECaps functions inside to create the black sheet in the top and bottom planes
	// Change the name: ConstructANNIECaps -> ConstructANNIECapsSheet to avoid conflicts with version in WCSimConstructANNIECylinder.cc

	ConstructANNIECapsSheet(-1);
	ConstructANNIECapsSheet(1);

	//--------------------------------------------------------------
	//--------------------ANNIE PMT holders ------------------------
	//--------------------------------------------------------------

	//If desired, construct ANNIE PMT holders
	//Geometrical properties of the holders are roughly extracted from the laser scan file (edges not clearly visible)

	G4bool HOLDER = WCSimTuningParams->GetHolder();
	G4cout <<"HOLDER variable: "<<HOLDER<<G4endl;
	G4double bsrff = WCSimTuningParams->GetBsrff();
	G4cout <<"Bsrff: "<<bsrff<<G4endl;

	if (HOLDER){
		ConstructANNIEHolders();
	}
	
	
	return logicWC;
	
}

//--------------------------------------------------------------------
//------------------------- Construct Caps ---------------------------
//--------------------------------------------------------------------

void WCSimDetectorConstruction::ConstructANNIECapsSheet(G4int zflip)
{
	//--------------------------------------------------------------------
	// ---------------------Add cap blacksheet----------------------------
	// -------------------------------------------------------------------
	
	G4double capBlackSheetZ[2] = {-WCBlackSheetThickness*zflip, 0.};
	G4double capBlackSheetRmin[2] = {0., 0.};
	G4double capBlackSheetRmax[2] = {WCIDRadius+WCBlackSheetThickness, 
									 WCIDRadius+WCBlackSheetThickness};
	G4VSolid* solidWCCapBlackSheet;
	if(WCBarrelRingNPhi*WCPMTperCellHorizontal == WCBarrelNumPMTHorizontal){
		solidWCCapBlackSheet
			= new G4Polyhedra("WCCapBlackSheet",  // name
							  0.*deg,             // phi start
							  totalAngle,         // total phi
							  WCBarrelRingNPhi,   // NPhi-gon
							  2,                  // z-planes
							  capBlackSheetZ,     // position of the Z planes
							  capBlackSheetRmin,  // min radius at the z planes
							  capBlackSheetRmax   // max radius at the Z planes
							  );
		// G4cout << *solidWCCapBlackSheet << G4endl;
	}
	
	G4LogicalVolume* logicWCCapBlackSheet =
		new G4LogicalVolume(solidWCCapBlackSheet,
							G4Material::GetMaterial("Blacksheet"),
							"WCCapBlackSheet",
							0,0,0);
	
	capAssemblyHeight = mainAnnulusHeight/2+1*mm+WCBlackSheetThickness;
	
	G4VPhysicalVolume* physiWCCapBlackSheet =
		new G4PVPlacement(0,
						  G4ThreeVector(0.,0.,capAssemblyHeight*zflip+InnerStructureCentreOffset),
						  logicWCCapBlackSheet,
						  "WCCapBlackSheet",
						  logicWCBarrel,
						  false,
						  0,
						  true);
	G4cout<<"constructed cap blacksheet at height "<<capAssemblyHeight*zflip+InnerStructureCentreOffset<<G4endl;
	WCCylInfo[(zflip>0)]=(capAssemblyHeight*zflip+InnerStructureCentreOffset)/10.;
	
	G4LogicalBorderSurface * WaterBSBottomCapSurface 
			= new G4LogicalBorderSurface("WaterBSCapPolySurface",
										 physiWCBarrel,physiWCCapBlackSheet,
										 OpWaterBSSurface);
	
	G4VisAttributes* WCCapBlackSheetVisAtt 
			= new G4VisAttributes(G4Colour(0.9,0.2,0.2));
	
	// used for OGLSX
	if (Vis_Choice == "OGLSX"){
		G4VisAttributes* WCCapBlackSheetVisAtt 
			= new G4VisAttributes(G4Colour(0.9,0.2,0.2));
		if(!debugMode)
			logicWCCapBlackSheet->SetVisAttributes(G4VisAttributes::Invisible);
		else
			logicWCCapBlackSheet->SetVisAttributes(WCCapBlackSheetVisAtt);
	}
	// used for RayTracer (makes the caps blacksheet yellow)
	if (Vis_Choice == "RayTracer"){
		
		G4VisAttributes* WCCapBlackSheetVisAtt 
		= new G4VisAttributes(G4Colour(1.0,1.0,0.0));
		
		if(!debugMode)
				//logicWCCapBlackSheet->SetVisAttributes(G4VisAttributes::Invisible);
				logicWCCapBlackSheet->SetVisAttributes(WCCapBlackSheetVisAtt);
				// Use first line if you want to make the blacksheet on the caps invisible to view through
			else
				logicWCCapBlackSheet->SetVisAttributes(WCCapBlackSheetVisAtt);
	}
	
}

void WCSimDetectorConstruction::ConstructANNIEHolders(){

	G4cout <<"Construct ANNIE Holders"<<G4endl;

	//Tried to get rough dimensions of the ANNIE holder from the laser scan file
	//Not really sure about the thickness, assume 2cm thickness for now (should not be super important)
	G4Box *ANNIEHolder_Box = new G4Box("ANNIEHolder_Box",10.5*cm,17.75*cm,1.*cm);
	G4Tubs *ANNIEHolder_Tube = new G4Tubs("ANNIEHolder_Tube",0.0*cm,6.0*cm,1*cm,0*deg,360*deg);

	//Create combined logical volume of the Box + Tube to get holder with hole (Subtraction Solid)

	G4SubtractionSolid *solidANNIEHolder = new G4SubtractionSolid("ANNIEHolder",
																ANNIEHolder_Box,
																ANNIEHolder_Tube,
																0,
																G4ThreeVector(0.,0.,0.));

	//Check the material of the ANNIE holders somewhere!
	G4LogicalVolume *logANNIEHolder = new G4LogicalVolume(solidANNIEHolder,
			G4Material::GetMaterial("PVC"),
			"WCANNIEHolder",
			0,0,0);

	G4double dist_pmt_holder = 10.84;		//Holder is 20cm away from the front face of the ANNIE PMTs, WCSim center is 9.16cm away from front --> 10.84cm distance

	//Read in PMT positions again, and project position of holder positions from the PMT positions

	// Create rotation matrices for the orientations of the PMTs
	std::vector<G4RotationMatrix*> holder_rotation_matrices;
	// Bottom PMTs have panel number 0
	G4RotationMatrix *WCBottomCapRotation = new G4RotationMatrix();
	holder_rotation_matrices.push_back(WCBottomCapRotation);
	G4RotationMatrix* WCPMTRotation = new G4RotationMatrix;
	WCPMTRotation->rotateY(90.*deg);	
	// Barrel PMTs have panel numbers 1-8
	for(int facei=0; facei<WCBarrelRingNPhi; facei++){
		G4RotationMatrix* WCPMTRotationNext = new G4RotationMatrix(*WCPMTRotation);
		WCPMTRotationNext->rotateX((dPhi*facei)-67.5*deg+180*deg);
		holder_rotation_matrices.push_back(WCPMTRotationNext);
	}
	// Top PMTs have panel number 9
	G4RotationMatrix *WCTopCapRotation = new G4RotationMatrix();
	WCTopCapRotation->rotateY(180.*deg);
	holder_rotation_matrices.push_back(WCTopCapRotation);

	//Select only ANNIE PMTs and propagate their position outwards to get central holder position
	std::ifstream pmt_position_file("PMTPositions_Scan.txt");
	std::string next_pmt;
	double pmt_x, pmt_y, pmt_z, pmt_dirx, pmt_diry, pmt_dirz;
	double holder_x, holder_y, holder_z;
	int panel_nr, pmt_type;
	int HolderID;
	while (!pmt_position_file.eof()){
		pmt_position_file >> HolderID >> panel_nr >> pmt_x >> pmt_y >> pmt_z >> pmt_dirx >> pmt_diry >> pmt_dirz >> pmt_type;
		if (pmt_position_file.eof()) break;
		G4cout << "Read in PMT "<<HolderID<<", panel nr: "<<panel_nr<<", Position ("<<pmt_x<<","<<pmt_y<<","<<pmt_z<<"), PMT type: "<<pmt_type<<G4endl;
		
		if (pmt_type == 2){	//select only ANNIE (Hamamatsu 8inch) PMTs for the holders
		G4RotationMatrix *holder_rot = holder_rotation_matrices.at(panel_nr);

		//Shift the PMT position outwards
		pmt_x -= (pmt_dirx*dist_pmt_holder);
		pmt_y -= (pmt_diry*dist_pmt_holder);
		pmt_z -= (pmt_dirz*dist_pmt_holder);

		holder_x = pmt_x*cm;
		holder_y = (168.1-pmt_z)*cm;
		holder_z = ((pmt_y+14.45))*cm;
		//pmt_z_shift = ((pmt_y+14.45)-InnerStructureCentreOffset/10.)*cm;
		G4cout <<"Edited Holder position ("<<holder_x<<","<<holder_y<<","<<holder_z<<")"<<G4endl;
		
		G4ThreeVector HolderPosition(holder_x,holder_y,holder_z);
		G4VPhysicalVolume *physicalHolder = new G4PVPlacement(holder_rot,	//its rotation
															HolderPosition,		//its position
															logANNIEHolder,			//its logical volume
															"ANNIE-Holder",			//its name
															logicWCBarrel,		//its mother volume
															false,				//no boolean operations
															HolderID,				//ID for this PMT (=channelkey in data)
															true);				//check overlaps*/
		G4LogicalBorderSurface* ANNIEHolderSurface
                = new G4LogicalBorderSurface("ANNIEHolderSurface",
                                                                         physiWCBarrel,
                                                                         physicalHolder,
                                                                         HolderOpSurface);
		}
	}
	pmt_position_file.close();

}
