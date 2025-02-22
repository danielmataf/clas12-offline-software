/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.jlab.geom.detector.alert.AHDC;
import org.jlab.geom.DetectorId;
import org.jlab.geom.abs.AbstractSector;
/**
 *
 * @author sergeyeva
 */
public class MYSector_ALERTDCWire extends AbstractSector<MYSuperlayer_ALERTDCWire> {
    protected MYSector_ALERTDCWire(int sectorId) {
        super(DetectorId.DC, sectorId);
    }
    
    /**
     * Returns "ALERT DC Sector".
     * @return "ALERT DC Sector"
     */
    @Override
    public String getType() {
        return "ALERT DC Sector";
    }
    
}
