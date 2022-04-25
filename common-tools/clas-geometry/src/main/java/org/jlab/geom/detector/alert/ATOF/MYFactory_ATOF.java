package org.jlab.geom.detector.alert.ATOF;

import java.util.ArrayList;
import java.util.List;
import org.jlab.geom.base.ConstantProvider;
import org.jlab.geom.base.DetectorTransformation;
import org.jlab.geom.base.Factory;
import org.jlab.geom.component.ScintillatorPaddle;
import org.jlab.geom.prim.Plane3D;
import org.jlab.geom.prim.Point3D;
import org.jlab.geom.prim.Shape3D;
import org.jlab.geom.prim.Transformation3D;
import org.jlab.geom.prim.Triangle3D;

/**
 * A Low Energy Recoil Tracker (ALERT)
 * This is to implement ALERT Time of Flight detector (ATOF)
 * There are 15 {@code MYSector}
 * and 10 {@code MYSuperlayer} per ATOF Detector.  There are 2 {@code MYLayer}s
 * and each layer contains 60 {@code ScintillatorPaddles} grouped by 4 into
 * 15 sectors. 
 * @author sergeyev
 * Implemented following the CNDDetector example
 */
public class MYFactory_ATOF implements Factory <MYDetector, MYSector, MYSuperlayer, MYLayer> {
    @Override
    public MYDetector createDetectorCLAS(ConstantProvider cp) {
        return createDetectorSector(cp);
    }
    
    @Override
    public MYDetector createDetectorSector(ConstantProvider cp) {
        return createDetectorTilted(cp);
    }
    
    @Override
    public MYDetector createDetectorTilted(ConstantProvider cp) {
        return createDetectorLocal(cp);
    }
    
    // From here and on, I take CND example to test how I can create scintillator bars 
    // Creating new objects: detector, sector, superlayer, layer from here
    // say how much sectors I want
    int nsectors = 15;
    @Override
    public MYDetector createDetectorLocal(ConstantProvider cp) {
        MYDetector detector = new MYDetector();
        for (int sectorId=0; sectorId<nsectors; sectorId++)
            detector.addSector(createSector(cp, sectorId));
        return detector;
    }

    // say how much superlayers I want in one sector
    int nsuperl = 10;
    @Override
    public MYSector createSector(ConstantProvider cp, int sectorId) {
        if(!(0<=sectorId && sectorId<nsectors))
            throw new IllegalArgumentException("Error: invalid sector="+sectorId);
        MYSector sector = new MYSector(sectorId);
        for (int superlayerId=0; superlayerId<nsuperl; superlayerId++)
            sector.addSuperlayer(createSuperlayer(cp, sectorId, superlayerId));
        return sector;
    }

    // say how much layers I want in one superlayer (so, in one sector)
    int nlayers = 2;
    @Override
    public MYSuperlayer createSuperlayer(ConstantProvider cp, int sectorId, int superlayerId) {
        if(!(0<=sectorId && sectorId<nsectors))
            throw new IllegalArgumentException("Error: invalid sector="+sectorId);
        if(!(0<=superlayerId && superlayerId<nsuperl))
            throw new IllegalArgumentException("Error: invalid superlayer="+superlayerId);
        MYSuperlayer superlayer = new MYSuperlayer(sectorId, superlayerId);
        for (int layerId=0; layerId<nlayers; layerId++)
            superlayer.addLayer(createLayer(cp, sectorId, superlayerId, layerId));
        return superlayer;
    }
    
    // here I create a layer, using a set of elementary objects that are called "component"; 
    // here, "component" = ScintillatorPaddle;
    @Override
    public MYLayer createLayer(ConstantProvider cp, int sectorId, int superlayerId, int layerId) {
        if(!(0<=sectorId && sectorId<nsectors))
            throw new IllegalArgumentException("Error: invalid sector="+sectorId);
        if(!(0<=superlayerId && superlayerId<nsuperl))
            throw new IllegalArgumentException("Error: invalid superlayer="+superlayerId);
        if(!(0<=layerId && layerId<nlayers))
            throw new IllegalArgumentException("Error: invalid layer="+layerId);

        double R0       = 75.0d; // InnerRadius layer 0 units? mm?
        double R1       = 78.0d; // InnerRadius layer 1
        double dR0       = 3.0d; // Thickness layer 0
        double dR1       = 20.0d;// Thickness layer 1
        double gp       = 0.0d; // LateralGap =gap paddle
        double gl       = 0.0d; // AzimuthalGap =gap layer
        double phi      = Math.toRadians(24.0d); // OpenAngle of a sector or block of paddles
        double len      = 30.0d; // Length in Z 
        double widthT   = 10.27d;// HigherBase of one paddle
        double widthM   = 8.17d;// MiddleBase of one paddle interface between layer 0 and 1
        double widthB   = 7.86d;// LowerBase of one paddle
        double z        = 0.0d; // UpstreamZOffset
        
        MYLayer layer = new MYLayer(sectorId, superlayerId, layerId);
        
        List<Plane3D> planes = new ArrayList();
    
        // we can have components grouped in blocks, here we have 4 component per 1 "block";
        
        int component = sectorId*4;
        double len_b = superlayerId*len;
        double len_f = (superlayerId+1)*len;
        double Rl = R0;
        double dR = dR0;
        double widthTl = widthM;
        double widthBl = widthB;
        
        double alpha   = Math.toRadians(3);
                
        if (layerId==1) {Rl = R1;
        dR = dR1; 
        widthTl = widthT;
        widthBl =widthM;
        }
        
        //for (int block=0; block<nblocks; block++) {
            // Left internal Paddle, with right edge at x = 0
            /*
            Point3D p0L = new Point3D(0,      0,  len_f);
            Point3D p1L = new Point3D(widthBl, 0,  len_f);
            Point3D p2L = new Point3D(widthTl, dR, len_f);
            Point3D p3L = new Point3D(0,      dR, len_f);
            Point3D p4L = new Point3D(0,      0,  len_b);
            Point3D p5L = new Point3D(widthBl, 0,  len_b);
            Point3D p6L = new Point3D(widthTl, dR, len_b);
            Point3D p7L = new Point3D(0,      dR, len_b);
            ScintillatorPaddle lPaddle = new ScintillatorPaddle(component, p0L, p1L, p2L, p3L, p4L, p5L, p6L, p7L);
            */
            Point3D p0L = new Point3D(0,      0,  len_f);
            Point3D p1L = new Point3D(widthBl*Math.cos(2*alpha), -widthBl*Math.sin(2*alpha),  len_f);
            Point3D p2L = new Point3D(widthTl*Math.cos(2*alpha), dR/Math.cos(alpha)-widthTl*Math.sin(2*alpha), len_f);
            Point3D p3L = new Point3D(0,      dR/Math.cos(alpha), len_f);
            
            Point3D p4L = new Point3D(0,      0,  len_b);
            Point3D p5L = new Point3D(widthBl*Math.cos(2*alpha), -widthBl*Math.sin(2*alpha),  len_b);
            Point3D p6L = new Point3D(widthTl*Math.cos(2*alpha), dR/Math.cos(alpha)-widthTl*Math.sin(2*alpha), len_b);
            Point3D p7L = new Point3D(0,      dR/Math.cos(alpha), len_b);
            ScintillatorPaddle lPaddle = new ScintillatorPaddle(component, p0L, p1L, p2L, p3L, p4L, p5L, p6L, p7L);
            
            //component = block==0? 0 : component+1;
            component = component+1;
            // Right internal Paddle, with left edge at x = 0
            /*
            Point3D p0R = new Point3D(-widthBl, 0,  len_f);
            Point3D p1R = new Point3D(0,       0,  len_f);
            Point3D p2R = new Point3D(0,       dR, len_f);
            Point3D p3R = new Point3D(-widthTl, dR, len_f);
            Point3D p4R = new Point3D(-widthBl, 0,  len_b);
            Point3D p5R = new Point3D(0,       0,  len_b);
            Point3D p6R = new Point3D(0,       dR, len_b);
            Point3D p7R = new Point3D(-widthTl, dR, len_b);
            ScintillatorPaddle rPaddle = new ScintillatorPaddle(component, p0R, p1R, p2R, p3R, p4R, p5R, p6R, p7R);
            */
            Point3D p0R = new Point3D(-widthBl*Math.cos(2*alpha), -widthBl*Math.sin(2*alpha),  len_f);
            Point3D p1R = new Point3D(0,       0,  len_f);
            Point3D p2R = new Point3D(0,      dR/Math.cos(alpha), len_f);
            Point3D p3R = new Point3D(-widthTl*Math.cos(2*alpha), dR/Math.cos(alpha)-widthTl*Math.sin(2*alpha), len_f);
            
            Point3D p4R = new Point3D(-widthBl*Math.cos(2*alpha), -widthBl*Math.sin(2*alpha),  len_b);
            Point3D p5R = new Point3D(0,       0,  len_b);
            Point3D p6R = new Point3D(0,      dR/Math.cos(alpha), len_b);
            Point3D p7R = new Point3D(-widthTl*Math.cos(2*alpha), dR/Math.cos(alpha)-widthTl*Math.sin(2*alpha), len_b);
            ScintillatorPaddle rPaddle = new ScintillatorPaddle(component, p0R, p1R, p2R, p3R, p4R, p5R, p6R, p7R);
            component = component+1;
            
            //here add external left & right paddles in order to have 4 paddles per sector
            /*
            Point3D p0Lext = new Point3D(widthBl,      0,  len_f);
            Point3D p1Lext = new Point3D(widthBl*2, 0,  len_f);
            Point3D p2Lext = new Point3D(widthTl*2, dR, len_f);
            Point3D p3Lext = new Point3D(widthTl, dR, len_f);
            Point3D p4Lext = new Point3D(widthBl,      0,  len_b);
            Point3D p5Lext = new Point3D(widthBl*2, 0,  len_b);
            Point3D p6Lext = new Point3D(widthTl*2, dR, len_b);
            Point3D p7Lext = new Point3D(widthTl,      dR, len_b);
            ScintillatorPaddle lextPaddle = new ScintillatorPaddle(component, p0Lext, p1Lext, p2Lext, p3Lext, p4Lext, p5Lext, p6Lext, p7Lext);
            */
            Point3D p0Lext = new Point3D(widthBl*Math.cos(2*alpha), -widthBl*Math.sin(2*alpha),  len_f);
            Point3D p1Lext = new Point3D(widthBl*Math.cos(2*alpha) + widthBl*Math.cos(4*alpha), -widthBl*Math.sin(2*alpha) - widthBl*Math.sin(4*alpha),  len_f);
            Point3D p2Lext = new Point3D(widthTl*Math.cos(2*alpha) + widthTl*Math.cos(4*alpha), dR/Math.cos(alpha)-widthTl*Math.sin(2*alpha) - widthTl*Math.sin(4*alpha), len_f);
            Point3D p3Lext = new Point3D(widthTl*Math.cos(2*alpha), dR/Math.cos(alpha)-widthTl*Math.sin(2*alpha), len_f);
            
            Point3D p4Lext = new Point3D(widthBl*Math.cos(2*alpha), -widthBl*Math.sin(2*alpha),  len_b);
            Point3D p5Lext = new Point3D(widthBl*Math.cos(2*alpha) + widthBl*Math.cos(4*alpha), -widthBl*Math.sin(2*alpha) - widthBl*Math.sin(4*alpha),  len_b);
            Point3D p6Lext = new Point3D(widthTl*Math.cos(2*alpha) + widthTl*Math.cos(4*alpha), dR/Math.cos(alpha)-widthTl*Math.sin(2*alpha) - widthTl*Math.sin(4*alpha), len_b);
            Point3D p7Lext = new Point3D(widthTl*Math.cos(2*alpha), dR/Math.cos(alpha)-widthTl*Math.sin(2*alpha), len_b);
            ScintillatorPaddle lextPaddle = new ScintillatorPaddle(component, p0Lext, p1Lext, p2Lext, p3Lext, p4Lext, p5Lext, p6Lext, p7Lext);
            component = component+1;
            
            /*
            Point3D p0Rext = new Point3D(-widthBl*2, 0,  len_f);
            Point3D p1Rext = new Point3D(-widthBl,       0,  len_f);
            Point3D p2Rext = new Point3D(-widthTl,       dR, len_f);
            Point3D p3Rext = new Point3D(-widthTl*2, dR, len_f);
            Point3D p4Rext = new Point3D(-widthBl*2, 0,  len_b);
            Point3D p5Rext = new Point3D(-widthBl,       0,  len_b);
            Point3D p6Rext = new Point3D(-widthTl,       dR, len_b);
            Point3D p7Rext = new Point3D(-widthTl*2, dR, len_b);
            ScintillatorPaddle rextPaddle = new ScintillatorPaddle(component, p0Rext, p1Rext, p2Rext, p3Rext, p4Rext, p5Rext, p6Rext, p7Rext);
            */
            Point3D p0Rext = new Point3D(-widthBl*Math.cos(2*alpha) - widthBl*Math.cos(4*alpha), -widthBl*Math.sin(2*alpha) - widthBl*Math.sin(4*alpha),  len_f);
            Point3D p1Rext = new Point3D(-widthBl*Math.cos(2*alpha), -widthBl*Math.sin(2*alpha),  len_f);
            Point3D p2Rext = new Point3D(-widthTl*Math.cos(2*alpha), dR/Math.cos(alpha)-widthTl*Math.sin(2*alpha), len_f);
            Point3D p3Rext = new Point3D(-widthTl*Math.cos(2*alpha) - widthTl*Math.cos(4*alpha), dR/Math.cos(alpha)-widthTl*Math.sin(2*alpha) - widthTl*Math.sin(4*alpha), len_f);
            
            Point3D p4Rext = new Point3D(-widthBl*Math.cos(2*alpha) - widthBl*Math.cos(4*alpha), -widthBl*Math.sin(2*alpha) - widthBl*Math.sin(4*alpha),  len_b);
            Point3D p5Rext = new Point3D(-widthBl*Math.cos(2*alpha), -widthBl*Math.sin(2*alpha),  len_b);
            Point3D p6Rext = new Point3D(-widthTl*Math.cos(2*alpha), dR/Math.cos(alpha)-widthTl*Math.sin(2*alpha), len_b);
            Point3D p7Rext = new Point3D(-widthTl*Math.cos(2*alpha) - widthTl*Math.cos(4*alpha), dR/Math.cos(alpha)-widthTl*Math.sin(2*alpha) - widthTl*Math.sin(4*alpha), len_b);
            ScintillatorPaddle rextPaddle = new ScintillatorPaddle(component, p0Rext, p1Rext, p2Rext, p3Rext, p4Rext, p5Rext, p6Rext, p7Rext);

            // Move the paddles into position relative to eachother
            lPaddle.translateXYZ( gp/2, 0, 0);
            rPaddle.translateXYZ(-gp/2, 0, 0);
            // adding the 3rd & 4th paddles to create a block of 4 paddles instead of 2 paddles
            lextPaddle.translateXYZ( gp/2, 0, 0);
            rextPaddle.translateXYZ(-gp/2, 0, 0);
            
            // Translate the paddles up to their proper radial distance and
            // along the z-axis
            //double r = Rl + (dR + gl) * layerId;
            double r = Rl;
            lPaddle.translateXYZ(0, r, -z);
            rPaddle.translateXYZ(0, r, -z);
            // adding the 3rd & 4th paddles to create a block of 4 paddles instead of 2 paddles
            lextPaddle.translateXYZ(0, r, -z);
            rextPaddle.translateXYZ(0, r, -z);
            
            // Rotate the paddles into their final position
            //double theta = block*phi - Math.toRadians(90);
            double theta = sectorId*phi - Math.toRadians(90);
            lPaddle.rotateZ(theta);
            rPaddle.rotateZ(theta);
            // adding the 3rd & 4th paddles to create a block of 4 paddles instead of 2 paddles
            lextPaddle.rotateZ(theta);
            rextPaddle.rotateZ(theta);
            
            // Add the paddles to the list
            layer.addComponent(lPaddle);
            layer.addComponent(rPaddle);
            layer.addComponent(lextPaddle);
            layer.addComponent(rextPaddle);
            
            Plane3D plane = new Plane3D(0, r, 0, 0, 1, 0);
            plane.rotateZ(theta);
            planes.add(plane);
        //}
        
        /*
        Plane3D fPlane = new Plane3D(0, 0, -z,     0, 0, 1);
        Plane3D bPlane = new Plane3D(0, 0, -z+len, 0, 0, 1);
        Shape3D boundary = layer.getBoundary();
        for (int block=0; block<nblocks; block++) {
            Plane3D plane0 = planes.get(block);
            Plane3D plane1 = planes.get((block+1)%nblocks);
            Plane3D plane2 = planes.get((block+2)%nblocks);
            Point3D p0 = new Point3D();
            Point3D p1 = new Point3D();
            Point3D p2 = new Point3D();
            Point3D p3 = new Point3D();
            fPlane.intersection(plane0, plane1, p0);
            bPlane.intersection(plane0, plane1, p1);
            fPlane.intersection(plane1, plane2, p2);
            bPlane.intersection(plane1, plane2, p3);
            Triangle3D f1 = new Triangle3D(p0, p1, p2);
            Triangle3D f2 = new Triangle3D(p3, p2, p1);
            boundary.addFace(f1);
            boundary.addFace(f2);
        }
        
        layer.getPlane().copy(fPlane);
        */
        return layer;
    }

    /**
     * Returns "CND Factory".
     * @return "CND Factory"
     */
    @Override
    public String getType() {
        return "CND Factory";
    }

    @Override
    public void show() {
        System.out.println(this);
    }
    
    @Override
    public String toString() {
        return getType();
    }

    @Override
    public Transformation3D getTransformation(ConstantProvider cp, int sector, int superlayer, int layer) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public DetectorTransformation getDetectorTransform(ConstantProvider cp) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

   
}

