package org.jlab.rec.ahdc.Banks;

import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.rec.ahdc.HelixFit.HelixFitObject;
import org.jlab.rec.ahdc.Hit.Hit;
import java.util.List;

public class RecoBankWriter {

  public DataBank fillAHDCHitsBank(DataEvent event, List<Hit> hitlist) {
    if (hitlist == null) {
      return null;
    }
    if (hitlist.size() == 0) {
      return null;
    }

    DataBank bank = event.createBank("AHDCRec::Hits", hitlist.size());

    for (int i = 0; i < hitlist.size(); i++) {

      bank.setShort("ID", i, (short) hitlist.get(i).get_Id());
      bank.setByte("layer", i, (byte) hitlist.get(i).get_Layer());
      bank.setByte("superlayer", i, (byte) hitlist.get(i).get_Super_layer());
      bank.setInt("wire", i, hitlist.get(i).get_Wire());
      bank.setDouble("Doca", i, hitlist.get(i).get_Doca());
    }

    return bank;
  }

  public DataBank fillAHDCMCTrackBank(DataEvent event) {

    DataBank particle = event.getBank("MC::Particle");
    double x_mc = particle.getFloat("vx", 0);
    double y_mc = particle.getFloat("vy", 0);
    double z_mc = particle.getFloat("vz", 0);
    double px_mc = particle.getFloat("px", 0) * 1000;
    double py_mc = particle.getFloat("py", 0) * 1000;
    double pz_mc = particle.getFloat("pz", 0) * 1000;

    int row = 0;
    DataBank bank = event.createBank("AHDCRec::MC", row + 1);
    bank.setFloat("x", row, (float) x_mc);
    bank.setFloat("y", row, (float) y_mc);
    bank.setFloat("z", row, (float) z_mc);
    bank.setFloat("px", row, (float) px_mc);
    bank.setFloat("py", row, (float) py_mc);
    bank.setFloat("pz", row, (float) pz_mc);

    return bank;
  }

  public DataBank fillAHDCTrackBank(DataEvent event, HelixFitObject ho) {

    double x = ho.get_X0();
    double y = ho.get_Y0();
    double z = ho.get_Z0();
    double px = ho.get_px();
    double py = ho.get_py();
    double pz = ho.get_pz();

    int row = 0;
    DataBank bank = event.createBank("AHDCRec::Track", row + 1);
    bank.setFloat("x", row, (float) x);
    bank.setFloat("y", row, (float) y);
    bank.setFloat("z", row, (float) z);
    bank.setFloat("px", row, (float) px);
    bank.setFloat("py", row, (float) py);
    bank.setFloat("pz", row, (float) pz);

    return bank;
  }
}
