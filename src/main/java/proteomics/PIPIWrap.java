package proteomics;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.Index.BuildIndex;
import proteomics.PTM.FindPTM;
import proteomics.Search.CalEValue;
import proteomics.Search.CalSubscores;
import proteomics.Search.CalXcorr;
import proteomics.Search.Search;
import proteomics.Segment.InferenceSegment;
import proteomics.TheoSeq.MassTool;
import proteomics.Types.*;

import java.text.DecimalFormat;
import java.util.*;
import java.util.concurrent.Callable;

public class PIPIWrap implements Callable<List<FinalResultEntry>> {

    private static final Logger logger = LoggerFactory.getLogger(PIPIWrap.class);

    private final BuildIndex buildIndexObj;
    private final MassTool massToolObj;
    private final InferenceSegment inference3SegmentObj;
    private final Map<Integer, SpectrumEntry> numSpectrumMap;
    private final NavigableMap<Float, List<Integer>> subMassNumMap;
    private final Map<String, TreeSet<Integer>> siteMass1000Map;
    private final float ms1Tolerance;
    private final int ms1ToleranceUnit;
    private final float ms2Tolerance;
    private final float minPtmMass;
    private final float maxPtmMass;
    private final int maxMs2Charge;
    private final int batchStartIdx;
    private final Map<Integer, String> pepTruth;
    private final Map<Integer, String> newTruth;

    public PIPIWrap(BuildIndex buildIndexObj, MassTool massToolObj, InferenceSegment inference3SegmentObj, Map<Integer, SpectrumEntry> numSpectrumMap, NavigableMap<Float, List<Integer>> subMassNumMap, Map<String, TreeSet<Integer>> siteMass1000Map, float ms1Tolerance, int ms1ToleranceUnit, float ms2Tolerance, float minPtmMass, float maxPtmMass, int maxMs2Charge, int batchStartIdx, Map<Integer, String> pepTruth, Map<Integer, String> newTruth) {
        this.buildIndexObj = buildIndexObj;
        this.massToolObj = massToolObj;
        this.inference3SegmentObj = inference3SegmentObj;
        this.numSpectrumMap = numSpectrumMap;
        this.subMassNumMap = subMassNumMap;
        this.siteMass1000Map = siteMass1000Map;
        this.ms1Tolerance = ms1Tolerance;
        this.ms1ToleranceUnit = ms1ToleranceUnit;
        this.ms2Tolerance = ms2Tolerance;
        this.minPtmMass = minPtmMass;
        this.maxPtmMass = maxPtmMass;
        this.maxMs2Charge = maxMs2Charge;
        this.batchStartIdx = batchStartIdx;
        this.pepTruth = pepTruth;
        this.newTruth = newTruth;
    }

    @Override
    public List<FinalResultEntry> call() {
        try {
            // Inference amino acids.
            Map<Integer, SparseVector> numCodeMap = new HashMap<>();
            Map<Integer, List<ThreeExpAA>> numExp3aaLists = new HashMap<>();
            for (List<Integer> numList : subMassNumMap.values()) {
                for (int scanNum : numList) {
                    SpectrumEntry spectrumEntry = numSpectrumMap.get(scanNum);

                    // Coding
                    List<ThreeExpAA> expAaLists = inference3SegmentObj.inferSegmentLocationFromSpectrum(scanNum, spectrumEntry);
                    if (spectrumEntry.scanNum == 1627){
                        System.out.println("here");
                    }
                    if (!expAaLists.isEmpty()) {
                        numExp3aaLists.put(scanNum, expAaLists);
                        numCodeMap.put(scanNum, inference3SegmentObj.generateSegmentIntensityVector(expAaLists));
                    }
                }
            }

            // Begin search.
            Search searchObj = new Search(buildIndexObj, numCodeMap, inference3SegmentObj, subMassNumMap, ms1Tolerance, ms1ToleranceUnit, minPtmMass, maxPtmMass, maxMs2Charge, batchStartIdx, pepTruth);
            Map<Integer, List<Peptide>> ptmOnlyMap = new HashMap<>();
            for (int scanNum : searchObj.getPTMOnlyResult().keySet()){
                List<Peptide> newPTMOnlyList = new LinkedList<>();
                for (Peptide oldPep : searchObj.getPTMOnlyResult().get(scanNum)){
                    newPTMOnlyList.add((Peptide)oldPep.clone());
                }
                ptmOnlyMap.put(scanNum, newPTMOnlyList);
            }
            Map<Integer, List<Peptide>> ptmFreeMap = new HashMap<>();
            for (int scanNum : searchObj.getPTMFreeResult().keySet()){
                List<Peptide> newPTMFreeList = new LinkedList<>();
                for (Peptide oldPep : searchObj.getPTMFreeResult().get(scanNum)){
                    newPTMFreeList.add((Peptide)oldPep.clone());
                }
                ptmFreeMap.put(scanNum, newPTMFreeList);
            }
            Map<Integer, Double> truthScore = searchObj.getTruthScore();
            DecimalFormat df = new DecimalFormat("0.000");
            int numRank = 0;
            int numRank1 = 0;
            int numCorrect = 0;
            for (int scanNum : pepTruth.keySet()) { // only check the accuracy for scans in simulation 1
                if (!truthScore.containsKey(scanNum)){
                    continue;
                }
                if (scanNum == 2537 || scanNum == 2095 || scanNum == 2052 || scanNum == 22598){
                    for (Peptide pep : ptmOnlyMap.get(scanNum)){
                        System.out.println("ptmOnly "+pep.getPTMFreeSeq()+": "+pep.getNormalizedCrossCorr());
                    }

                    for (Peptide pep : ptmFreeMap.get(scanNum)){
                        System.out.println("ptmFree "+pep.getPTMFreeSeq()+": "+pep.getNormalizedCrossCorr());
                    }
                }

                TreeMap<Float, String> scanRes = new TreeMap<>();
                if (ptmOnlyMap.containsKey(scanNum)) {
                    for (Peptide candidate : ptmOnlyMap.get(scanNum)) {
                        scanRes.put((float) candidate.getNormalizedCrossCorr(), candidate.getPTMFreeSeq());
                    }
                }
                if (ptmFreeMap.containsKey(scanNum)) {
                    for (Peptide candidate : ptmFreeMap.get(scanNum)) {
                        scanRes.put((float) candidate.getNormalizedCrossCorr(), candidate.getPTMFreeSeq());
                    }
                }

                String toPrint = scanNum + " " + pepTruth.get(scanNum) + ": ";
                int truthRank = -1;
                int i = 0;
                for (Map.Entry<Float, String> candi : scanRes.entrySet()){
                    if (candi.getValue().equals(pepTruth.get(scanNum))) {
                        truthRank = scanRes.size()-i;
                    }
                    if (i >= scanRes.size()-3){
                        toPrint += candi.getValue() +"="+candi.getKey()+"  ";
                    }
                    i++;
                }
                numRank++;
                double topScore = 0;
                String topSeq = "";
                for (double s : scanRes.keySet()){
                    if (s>=topScore){
                        topScore =s ;
                        topSeq = scanRes.get((float) s);
                    }
                }
                if (truthRank == 1){
//                    System.out.println(scanNum + " ,top1, "+ topScore);
                    numRank1++;
                }
//                System.out.println(scanNum + " , '"+ topSeq+  "' , "+ topSeq.equals(pepTruth.get(scanNum)));
            }
//            System.out.println(numRank1+" "+numRank);
            Map<Integer, List<Peptide>> numCandidateMapNoPTMStep = new HashMap<>(ptmFreeMap);
            for (int scanNum : ptmOnlyMap.keySet()) {
                if (scanNum == 13701 || scanNum == 2095 || scanNum == 2052 || scanNum == 22598){
                    System.out.println("here");
                }
                if (numCandidateMapNoPTMStep.containsKey(scanNum)) {
                    numCandidateMapNoPTMStep.get(scanNum).addAll(ptmOnlyMap.get(scanNum));
                } else {
                    numCandidateMapNoPTMStep.put(scanNum, ptmOnlyMap.get(scanNum));
                }
            }
            logger.debug("Analyzing PTMs...");
            FindPTM findPtmObj = new FindPTM(searchObj.getPTMOnlyResult(), numSpectrumMap, numExp3aaLists, massToolObj, siteMass1000Map, minPtmMass, maxPtmMass, ms1Tolerance, ms1ToleranceUnit, ms2Tolerance, batchStartIdx);
            Map<Integer, List<Peptide>> ptmOnlyTemp = findPtmObj.getPeptidesWithPTMs();

            logger.debug("Calculating final score...");
            Map<Integer, List<Peptide>> numCandidateMap = new HashMap<>(searchObj.getPTMFreeResult());
            for (int scanNum : ptmOnlyTemp.keySet()) {
                if (scanNum == 13701 || scanNum == 2095 || scanNum == 2052 || scanNum == 22598){
                    System.out.println("here");
                }
                if (numCandidateMap.containsKey(scanNum)) {
                    numCandidateMap.get(scanNum).addAll(ptmOnlyTemp.get(scanNum));
                } else {
                    numCandidateMap.put(scanNum, ptmOnlyTemp.get(scanNum));
                }
            }

            CalXcorr calXcorrObj = new CalXcorr(numCandidateMap, numSpectrumMap, massToolObj, buildIndexObj, numCandidateMapNoPTMStep);
            List<FinalResultEntry> subScoredPsms = calXcorrObj.getScoredPSMs();

            new CalSubscores(subScoredPsms, numSpectrumMap, ms2Tolerance);

            logger.debug("Estimating E-value for each PSM...");
            for (FinalResultEntry psm : subScoredPsms) {
                new CalEValue(psm);
            }
            return subScoredPsms;
        } catch (Exception ex) {
            logger.error(ex.getMessage());
            ex.printStackTrace();
            System.exit(1);
        }
        return new LinkedList<>();
    }
}
