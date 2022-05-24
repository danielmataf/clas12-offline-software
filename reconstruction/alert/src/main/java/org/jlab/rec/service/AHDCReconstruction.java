package org.jlab.rec.service;

import org.jlab.clas.reco.ReconstructionEngine;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataSource;
import org.jlab.io.hipo.HipoDataSync;
import org.jlab.rec.ahdc.Banks.RecoBankWriter;
import org.jlab.rec.ahdc.Cluster.Cluster;
import org.jlab.rec.ahdc.Cluster.ClusterFinder;
import org.jlab.rec.ahdc.Distance.Distance;
import org.jlab.rec.ahdc.HelixFit.HelixFitJava;
import org.jlab.rec.ahdc.HelixFit.HelixFitObject;
import org.jlab.rec.ahdc.Hit.Hit;
import org.jlab.rec.ahdc.Hit.HitReader;
import org.jlab.rec.ahdc.HoughTransform.HoughTransform;
import org.jlab.rec.ahdc.PreCluster.PreCluster;
import org.jlab.rec.ahdc.PreCluster.PreClusterFinder;
import org.jlab.rec.ahdc.Track.Track;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

public class AHDCReconstruction extends ReconstructionEngine {

  public AHDCReconstruction() {
    super("ALERT", "ouillon", "0.0.1");
  }

  @Override
  public boolean init() {
    return true;
  }

  /**
   * Command to run the engine : ./coatjava/bin/recon-util org.jlab.rec.service.AHDCReconstruction
   * org.jlab.clas.swimtools.MagFieldsEngine -i input.hipo -o output.hipo
   *
   * @param event Hipo DataEvent
   * @return 0
   */
  @Override
  public boolean processDataEvent(DataEvent event) {

    String method = "distance";

    if (event.hasBank("AHDC::adc")) {

      // VI) Write bank -------------------------------------------------------
      RecoBankWriter writer = new RecoBankWriter();
      DataBank recoBank = writer.fillAHDCTrackBank(event);

      event.appendBank(recoBank);
      // ----------------------------------------------------------------------

      // I) Read hit
      HitReader hitRead = new HitReader();
      hitRead.fetch_AHDCHits(event);
      List<Hit> AHDC_Hits = hitRead.get_AHDCHits();

      // II) Create PreCluster
      PreClusterFinder preclusterfinder = new PreClusterFinder();
      preclusterfinder.findPreCluster(AHDC_Hits);
      List<PreCluster> AHDC_PreClusters = preclusterfinder.get_AHDCPreClusters();

      // III) Create Cluster
      ClusterFinder clusterfinder = new ClusterFinder();
      clusterfinder.findCluster(AHDC_PreClusters);
      List<Cluster> AHDC_Clusters = clusterfinder.get_AHDCClusters();

      // IV) Track Finder
      List<Track> AHDC_Tracks = new ArrayList<>();
      if (method.equals("distance")) {
        // IV) a) Distance method
        Distance distance = new Distance();
        distance.find_track(AHDC_Clusters);
        AHDC_Tracks = distance.get_AHDCTracks();
      } else if (method.equals("hough")) {
        // IV) b) Hough Transform method
        HoughTransform houghtransform = new HoughTransform();
        houghtransform.find_tracks(AHDC_Clusters);
        AHDC_Tracks = houghtransform.get_AHDCTracks();
      }

      // V) Global fit
      for (Track track : AHDC_Tracks) {
        // Nombre de points
        int nbofpoint = 0;
        for (Cluster cluster : track.get_Clusters()) {
          nbofpoint++;
        }

        double[][] szPos = new double[nbofpoint][3];
        int j = 0;
        for (Cluster cluster : track.get_Clusters()) {
          szPos[j][0] = cluster.get_X();
          szPos[j][1] = cluster.get_Y();
          szPos[j][2] = cluster.get_Z();
          j++;
        }

        HelixFitJava h = new HelixFitJava();
        HelixFitObject ho = h.HelixFit(nbofpoint, szPos, 1);
        // Kalman Filter
        // Creation d'une helice pour le filtre de Kalman

        // Create a point on the beamline
        Cluster clus = new Cluster(ho.get_X0(), ho.get_Y0(), ho.get_Z0());
        track.get_Clusters().add(0, clus);

        // Nombre de points
        nbofpoint++;

        // szPos[nbofpoint][3] x y z pour chaque point
        szPos = new double[nbofpoint][3];
        j = 0;
        for (Cluster cluster : track.get_Clusters()) {
          szPos[j][0] = cluster.get_X();
          szPos[j][1] = cluster.get_Y();
          szPos[j][2] = cluster.get_Z();
          j++;
        }

        HelixFitJava h1 = new HelixFitJava();
        HelixFitObject ho1 = h1.HelixFit(nbofpoint, szPos, 0);


      }
    }
    return true;
  }

  public static void main(String[] args) {

    double starttime = System.nanoTime();

    String inputFile = "proton_10_100MeV.hipo";
    String outputFile = "output.hipo";

    if (new File(outputFile).delete()) System.out.println("output.hipo is delete.");

    System.err.println(" \n[PROCESSING FILE] : " + inputFile);

    AHDCReconstruction en = new AHDCReconstruction();
    en.init();

    HipoDataSource reader = new HipoDataSource();
    HipoDataSync writer = new HipoDataSync();
    reader.open(inputFile);
    writer.open(outputFile);
    while (reader.hasEvent()) {
      DataEvent event = reader.getNextEvent();
      en.processDataEvent(event);
      writer.writeEvent(event);
    }
    writer.close();

    System.out.println("finished " + (System.nanoTime() - starttime) * Math.pow(10, -9));
  }
}
