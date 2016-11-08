#ifndef __CONSTRUCT_MRD_PMT_VERBOSE__
#define __CONSTRUCT_MRD_PMT_VERBOSE__ 1
#endif

#include "WCSimDetectorConstruction.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Torus.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4Material.hh"
#include "G4Polycone.hh"
#include "G4PVPlacement.hh"
#include "G4LogicalBorderSurface.hh"

#include "G4SDManager.hh"
#include "WCSimWCSD.hh"
#include "WCSimPMTObject.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

//Flat faced PMT logical volume construction. Based on tank construction, but with different PMT geometry.

G4LogicalVolume* WCSimDetectorConstruction::ConstructFlatFacedPMT(G4String PMTName, G4String CollectionName, G4String detectorElement)
{
#ifdef __CONSTRUCT_MRD_PMT_VERBOSE__
  G4cout<<"Making PMTKey_t (PMTName="<<PMTName<<", CollectionName="<<CollectionName<<") "<<G4endl;
#endif
  PMTKey_t key(PMTName,CollectionName);
  PMTMap_t::iterator it = PMTLogicalVolumes.find(key);
  if (it != PMTLogicalVolumes.end()) {
#ifdef __CONSTRUCT_MRD_PMT_VERBOSE__
      G4cout<<"This key already exists in the PMTLogicalVolumes map. Restoring it."<<G4endl;
#endif
      //G4cout << "Restore PMT" << G4endl;
      return it->second;
  }
  //G4cout << "Create PMT" << G4endl;

#ifdef __CONSTRUCT_MRD_PMT_VERBOSE__
G4cout<<"Key not found; creating a new PMT logical volume"<<G4endl;
#endif

if (Vis_Choice == "RayTracer"){
    // Blue wireframe visual style
    // Used in the RayTracer visualizer
  G4VisAttributes* WCPMTVisAtt = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
  WCPMTVisAtt->SetForceSolid(true); // force the object to be visualized with a surface
  WCPMTVisAtt->SetForceAuxEdgeVisible(true); // force auxiliary edges to be shown 
}

else
   { // Gray wireframe visual style
    // used in OGLSX visualizer
  G4VisAttributes* WCPMTVisAtt = new G4VisAttributes(G4Colour(0.2,0.2,0.2));
  WCPMTVisAtt->SetForceWireframe(true);}

  G4double expose;
  G4double radius;
  G4double glassThickness;
#ifdef __CONSTRUCT_MRD_PMT_VERBOSE__
  G4cout<<"Getting pointer to PMT object corresponding to CollectionName "<<CollectionName<<G4endl;
#endif
  WCSimPMTObject *PMT = GetPMTPointer(CollectionName);
#ifdef __CONSTRUCT_MRD_PMT_VERBOSE__
  G4cout<<"Retrieved PMT with name "<<PMT->GetPMTName()<<G4endl;
#endif
  expose = PMT->GetExposeHeight();
  radius = PMT->GetRadius();
  glassThickness = PMT->GetPMTGlassThickness();

G4double shamferrad=radius/10.;	// arbitrary choice for now
G4double gelthickness=0.2*cm;		// equally arbitrary

#ifdef __CONSTRUCT_MRD_PMT_VERBOSE__
G4cout<<"Making the geometry"<<G4endl;
#endif

// start with a box of silicon 'optical grease' gel - this will smoosh against the light guide
G4Box* PMTgelbox = new G4Box("WCPMT",
                            radius+2.*cm,
                            radius+2.*cm,
                            (shamferrad+gelthickness+0.01*cm)/2.);

G4LogicalVolume* logicWCPMT = new G4LogicalVolume(PMTgelbox,
                            G4Material::GetMaterial("Silicone"),
                            "WCPMT",
                            0,0,0);

if (Vis_Choice == "RayTracer"){
// Makes the volume containing the PMT visible, solid, and forces the auxiliary edges to be viewed.
  G4VisAttributes* WCPMTVisAtt = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
  WCPMTVisAtt->SetForceSolid(true); // force the object to be visualized with a surface
  WCPMTVisAtt->SetForceAuxEdgeVisible(true); // force auxiliary edges to be shown 

    logicWCPMT->SetVisAttributes(WCPMTVisAtt);}

else{
// Makes the volume containg the PMT invisible for normal visualization
    logicWCPMT->SetVisAttributes(G4VisAttributes::Invisible);}

// now make up the PMT glass face from a series of parts
// first a torus to make up the chamfered edge
G4Torus* tmptorus = new G4Torus("tmptorus",
                            shamferrad-glassThickness,
                            shamferrad,
                            radius-shamferrad,
                            0,
                            2.0*pi);

// a solid disk to make up the main face of the PMT
G4Tubs* PMTface = new G4Tubs("PMTface",
                          0,
                          radius-shamferrad+0.2*cm,
                          (glassThickness/2.),
                          0,
                          2.0*pi);

// a solid disk to subtract the unwanted lower half of the torus
// note that subtraction solids should avoid aligned edges/faces, as precision loss could produce undefined results
G4Tubs* tmpsubtorusbottom = new G4Tubs("tmpsubtorusbottom",
                          0,
                          radius+1.*cm,
                          shamferrad,
                          0,
                          2.0*pi);

// one more solid disk to subtract the remaining unwanted inner quarter of the torus
G4Tubs* tmpsubtorusinner = new G4Tubs("tmpsubtorusinner",
                          0,
                          radius-shamferrad,
                          shamferrad+1.*cm,
                          0,
                          2.0*pi);

// now build up the PMT face from the components
G4ThreeVector zTrans = G4ThreeVector(0, 0, -shamferrad);
//G4RotationMatrix Rot = 0;
G4Transform3D transform = G4Transform3D(G4RotationMatrix(), zTrans);

// subtract the lower half of the torus
G4SubtractionSolid* tmpsub1 = new G4SubtractionSolid("tmpsub1",
                    tmptorus,
                    tmpsubtorusbottom,
                    transform);

// subtract the inner quarter of the torus
G4SubtractionSolid* PMTedge = new G4SubtractionSolid("PMTedge",
                    tmpsub1,
                    tmpsubtorusinner);

transform = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0,0,shamferrad-(glassThickness/2.)+0.02*cm));

// generate a union between the disk face and the remaining quarter torus acting as the chamfered edge
G4UnionSolid* roundedPMTface = new G4UnionSolid(CollectionName, 
                    PMTedge,
                    PMTface,
                    transform);

// make the logical volume to turn this PMT glass face into glass...
G4LogicalVolume* logicGlassFaceWCPMT = new G4LogicalVolume(roundedPMTface,
                    G4Material::GetMaterial("Glass"),
                    CollectionName,
                    0,0,0);

// ... and give it suitable visualization properties
G4VisAttributes* WCPMTVisAtt;
if (1){//Vis_Choice == "RayTracer"){
	WCPMTVisAtt = new G4VisAttributes(G4Colour(0,0.5,1.));
	WCPMTVisAtt->SetForceSolid(true); 					// force the object to be visualized with a surface
	WCPMTVisAtt->SetForceAuxEdgeVisible(true); 	// force auxiliary edges to be shown
} else { // Gray wireframe visual style used in OGLSX visualizer
	WCPMTVisAtt = new G4VisAttributes(G4Colour(0.2,0.2,0.2));
	WCPMTVisAtt->SetForceWireframe(true);
}


// For either visualization type, logicGlassFaceWCPMT will either be visible or invisible depending on which
// line is commented at the end of the respective if statements

//  if (Vis_Choice == "OGLSX")
//   { // Gray wireframe visual style
//    // used in OGLSX visualizer
//  G4VisAttributes* WCPMTVisAtt = new G4VisAttributes(G4Colour(0.2,0.2,0.2));
//  WCPMTVisAtt->SetForceWireframe(true);
//  //logicGlassFaceWCPMT->SetVisAttributes(G4VisAttributes::Invisible);
//  logicGlassFaceWCPMT->SetVisAttributes(WCPMTVisAtt);}

//  if (Vis_Choice == "RayTracer"){
//    // Blue wireframe visual style
//    // Used in the RayTracer visualizer
//  G4VisAttributes* WCPMTVisAtt = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
//  WCPMTVisAtt->SetForceSolid(true); // force the object to be visualized with a surface
//  WCPMTVisAtt->SetForceAuxEdgeVisible(true); // force auxiliary edges to be shown 
//  //logicGlassFaceWCPMT->SetVisAttributes(G4VisAttributes::Invisible);
//  logicGlassFaceWCPMT->SetVisAttributes(WCPMTVisAtt);}

// place this glass face into the box of optical grease
G4VPhysicalVolume* physiGlassFaceWCPMT = new G4PVPlacement(0,
                        G4ThreeVector(0, 0, -(shamferrad+gelthickness+0.01*cm)/2.),
                        logicGlassFaceWCPMT,
                        CollectionName,
                        logicWCPMT,
                        false,
                        0,
                        false);

/* ###### note: the solid, logical and physical volumes corresponding to the glass PMT face are all 
given the name CollectionName, which is passed in when ConstructPMT is called! ########### */

// make 'air' to go behind the glass cathode surface, as per other PMTs... dunno why it's not a vacuum.
G4Torus* tmpairtorus = new G4Torus("tmpairtorus",
                            0,
                            shamferrad-glassThickness,
                            radius-shamferrad,
                            0,
                            2.0*pi);

// subtract the bottom half of this torus
G4Tubs* tmpdisk = new G4Tubs("tmpdisk",
                          0,
                          radius+1.*cm,
                          shamferrad/2.,
                          0,
                          2.0*pi);

transform = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0,0,(-shamferrad/2.)));

G4SubtractionSolid* airtorus = new G4SubtractionSolid("airtorus",
                          tmpairtorus,
                          tmpdisk,
                          transform);

G4Tubs* airdisk = new G4Tubs("airdisk",
                          0,
                          radius-shamferrad,
                          shamferrad-glassThickness+0.005*cm,
                          0,
                          2.0*pi);

G4LogicalVolume* airtorus_log = new G4LogicalVolume(airtorus,
                    G4Material::GetMaterial("Air"),
                    "airtorus_log",
                    0,0,0);

G4LogicalVolume* airdisk_log = new G4LogicalVolume(airdisk,
                    G4Material::GetMaterial("Air"),
                    "airdisk_log",
                    0,0,0);

if (Vis_Choice == "RayTracer"){
// Adding color and forcing the inner portion of the PMT's to be viewed
  G4VisAttributes* WCPMTVisAtt = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
  WCPMTVisAtt->SetForceSolid(true); // force the object to be visualized with a surface
  WCPMTVisAtt->SetForceAuxEdgeVisible(true); // force auxiliary edges to be shown 

  airtorus_log->SetVisAttributes(WCPMTVisAtt);
  airdisk_log->SetVisAttributes(WCPMTVisAtt);}

else {
// Making the inner portion of the detector invisible for OGLSX visualization
  airtorus_log->SetVisAttributes(G4VisAttributes::Invisible);
  airdisk_log->SetVisAttributes(G4VisAttributes::Invisible);}

// place the air torus into the box of optical grease
G4VPhysicalVolume* airtorus_phys = new G4PVPlacement(0,
                        G4ThreeVector(0, 0, -(shamferrad+gelthickness+0.01*cm)/2.),
                        airtorus_log,
                        "airtorus_phys",
                        logicWCPMT,
                        false,
                        0,
                        false);

// place the air disk into the box of optical grease
G4VPhysicalVolume* airdisk_phys = new G4PVPlacement(0,
                        G4ThreeVector(0, 0, -(shamferrad+gelthickness+0.01*cm)/2.),
                        airdisk_log,
                        "airdisk_phys",
                        logicWCPMT,
                        false,
                        0,
                        false);

#ifdef __CONSTRUCT_MRD_PMT_VERBOSE__
  G4cout<<"Finished making the geometry"<<G4endl;
#endif
// end of physical construction
// ============================

  // Instantiate a new sensitive detector 
  // and register this sensitive detector volume with the SD Manager. 
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4String SDName = "/WCSim/";
  SDName += CollectionName;

#ifdef __CONSTRUCT_MRD_PMT_VERBOSE__
  G4cout<<"Searching for sensitive detector "<<SDName<<G4endl;
#endif
  // If there is no such sensitive detector with that SDName yet,
  // make a new one
  if( ! SDman->FindSensitiveDetector(SDName, false) ) {
#ifdef __CONSTRUCT_MRD_PMT_VERBOSE__
    G4cout<<"Sensitive detector not found. Making aWCPMT = new WCSimWCSD("<<CollectionName<<", "<<SDName
          <<", "<<"{DetectorConstruction}"<<", "<<detectorElement<<"), and adding to the SD manager"<<G4endl;
#endif
    aMRDPMT = new WCSimWCSD(CollectionName,SDName,this, detectorElement );
    SDman->AddNewDetector( aMRDPMT );
  }

  logicGlassFaceWCPMT->SetSensitiveDetector( aMRDPMT );

  PMTLogicalVolumes[key] = logicWCPMT;

#ifdef __CONSTRUCT_MRD_PMT_VERBOSE__
  G4cout<<"Adding optical surfaces to cathode"<<G4endl;
#endif
  //Add Logical Border Surface
  new G4LogicalBorderSurface("GlassCathodeSurface",
                             physiGlassFaceWCPMT,
                             airdisk_phys,
                             OpGlassCathodeSurface);

  //Add Logical Border Surface
  new G4LogicalBorderSurface("GlassCathodeSurface",
                             physiGlassFaceWCPMT,
                             airtorus_phys,
                             OpGlassCathodeSurface);

#ifdef __CONSTRUCT_MRD_PMT_VERBOSE__
  G4cout<<"returning logical volume"<<G4endl;
#endif
  return logicWCPMT;
}
