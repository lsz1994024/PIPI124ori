package proteomics.Search;


import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.Segment.InferenceSegment;
import proteomics.Index.BuildIndex;
import proteomics.TheoSeq.MassTool;
import proteomics.PIPI;
import proteomics.Types.Peptide;
import proteomics.Types.ResultEntry;
import proteomics.Types.SparseBooleanVector;
import proteomics.Types.SparseVector;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

public class Search {

    private static final Logger logger = LoggerFactory.getLogger(Search.class);
    private static final int rankNum = 10;

    private Map<Integer, List<Peptide>> ptmOnlyResult = new HashMap<>();
    private Map<Integer, List<Peptide>> ptmFreeResult = new HashMap<>();
    private final MassTool massToolObj;
    private final int maxMs2Charge;
    private Map<Integer, Double> truthScore = new HashMap<>();
    private Set<Integer> scanSearched = new HashSet<>();

    public Search(BuildIndex buildIndexObj, Map<Integer, SparseVector> numCodeMap, InferenceSegment inference3SegmentObj, NavigableMap<Float, List<Integer>> massNumMap, float ms1Tolerance, int ms1ToleranceUnit, float minPtmMass, float maxPtmMass, int maxMs2Charge, int batchStartIdx, Map<Integer, String> pepTruth) {
        massToolObj = buildIndexObj.returnMassToolObj();
        this.maxMs2Charge = maxMs2Charge;

        Map<String, Float> massTable = massToolObj.returnMassTable();
        Set<String> targetPeptideSet = buildIndexObj.returnPepMassMap().keySet();
        Map<String, Float> localPeptideMassMap = new HashMap<>();
        localPeptideMassMap.putAll(buildIndexObj.returnPepMassMap());
        localPeptideMassMap.putAll(buildIndexObj.returnDecoyPepMassMap());
        Map<String, Set<String>> peptideProteinMap = buildIndexObj.returnPepProMap();
        Map<String, String> proteinSeqMap = buildIndexObj.returnProPepMap();

        Map<Integer, Double> numCodeNormSquareMap = new HashMap<>();
        for (int scanNum : numCodeMap.keySet()) {
            SparseVector code = numCodeMap.get(scanNum);
            double threeNorm = code.norm2square();
            numCodeNormSquareMap.put(scanNum, threeNorm);
        }

        Map<Integer, PriorityQueue<ResultEntry>> ptmFreeResultMap = new HashMap<>();
        Map<Integer, PriorityQueue<ResultEntry>> ptmOnlyResultMap = new HashMap<>();

        for (String peptide : localPeptideMassMap.keySet()) {
            float peptideMass = localPeptideMassMap.get(peptide);
            SparseBooleanVector peptideThreeCode = inference3SegmentObj.generateSegmentBooleanVector(inference3SegmentObj.cutTheoSegment(peptide));
            double peptideCodeThreeNormSquare = peptideThreeCode.norm2square();

            float tolerance = ms1Tolerance;
            if (ms1ToleranceUnit == 1) {
                tolerance = peptideMass * ms1Tolerance / 1e6f;
            }

            float leftMass = Math.max(peptideMass + minPtmMass, massNumMap.firstKey());
            float rightMass = Math.min(peptideMass + maxPtmMass, massNumMap.lastKey());

            // The peptide mass is too small or too large.
            if (leftMass >= rightMass) {
                continue;
            }

            NavigableMap<Float, List<Integer>> subNumMap = massNumMap.subMap(leftMass, true, rightMass, true);

            if (subNumMap.isEmpty()) {
                continue;
            }

//            System.out.println(subNumMap.size());
            boolean isTarget = false;
            if (targetPeptideSet.contains(peptide)) {
                isTarget = true;
            }
            for (float scanMass : subNumMap.keySet()) {
                for (int scanNum : subNumMap.get(scanMass)) {
                    if (numCodeMap.containsKey(scanNum)) {
                        SparseVector scanCode = numCodeMap.get(scanNum);

                        if (scanNum == 1820 && peptide.equals("FTHAETHSSIK")){
                            System.out.println("here");
                        }
                        double scanNormSquare = numCodeNormSquareMap.get(scanNum);
                        double score = 0;
                        scanSearched.add(scanNum);
                        double temp1 = Math.sqrt(peptideCodeThreeNormSquare * scanNormSquare);
                        if (temp1 > 1e-6) {
                            score = peptideThreeCode.dot(scanCode) / temp1;
                            if (peptide.equals(pepTruth.get(scanNum))){
                                truthScore.put(scanNum, score);
//                                System.out.println(scanNum+ " "+ peptide + " "+score);
                            }
                        }
                        float deltaMass = scanMass - peptideMass;
                        if (isTarget) {
                            if (Math.abs(deltaMass) <= tolerance) {
                                // PTM-free
                                if (ptmFreeResultMap.containsKey(scanNum)) {
                                    PriorityQueue<ResultEntry> temp = ptmFreeResultMap.get(scanNum);
                                    if (temp.size() < rankNum) {
                                        temp.add(new ResultEntry(score, peptide, false, false));
                                    } else {
                                        if (score > temp.peek().score) {
                                            temp.poll();
                                            temp.add(new ResultEntry(score, peptide, false, false));
                                        }
                                    }
                                } else {
                                    PriorityQueue<ResultEntry> temp = new PriorityQueue<>(rankNum * 2);
                                    temp.add(new ResultEntry(score, peptide, false, false));
                                    ptmFreeResultMap.put(scanNum, temp);
                                }
                            }

                            if (Math.abs(deltaMass) > tolerance) {
                                // PTM-only
                                if (ptmOnlyResultMap.containsKey(scanNum)) {
                                    PriorityQueue<ResultEntry> temp = ptmOnlyResultMap.get(scanNum);
                                    if (temp.size() < rankNum) {
                                        temp.add(new ResultEntry(score, peptide, false, true));
                                    } else {
                                        if (score > temp.peek().score) {
                                            temp.poll();
                                            temp.add(new ResultEntry(score, peptide, false, true));
                                        }
                                    }
                                } else {
                                    PriorityQueue<ResultEntry> temp = new PriorityQueue<>(rankNum * 2);
                                    temp.add(new ResultEntry(score, peptide, false, true));
                                    ptmOnlyResultMap.put(scanNum, temp);
                                }
                            }
                        } else {
                            if (Math.abs(deltaMass) <= tolerance) {
                                // PTM-free
                                if (ptmFreeResultMap.containsKey(scanNum)) {
                                    PriorityQueue<ResultEntry> temp = ptmFreeResultMap.get(scanNum);
                                    if (temp.size() < rankNum) {
                                        temp.add(new ResultEntry(score, peptide, true, false));
                                    } else {
                                        if (score > temp.peek().score) {
                                            temp.poll();
                                            temp.add(new ResultEntry(score, peptide, true, false));
                                        }
                                    }
                                } else {
                                    PriorityQueue<ResultEntry> temp = new PriorityQueue<>(rankNum * 2);
                                    temp.add(new ResultEntry(score, peptide, true, false));
                                    ptmFreeResultMap.put(scanNum, temp);
                                }
                            }

                            if (Math.abs(deltaMass) > tolerance) {
                                // PTM-only
                                if (ptmOnlyResultMap.containsKey(scanNum)) {
                                    PriorityQueue<ResultEntry> temp = ptmOnlyResultMap.get(scanNum);
                                    if (temp.size() < rankNum) {
                                        temp.add(new ResultEntry(score, peptide, true, true));
                                    } else {
                                        if (score > temp.peek().score) {
                                            temp.poll();
                                            temp.add(new ResultEntry(score, peptide, true, true));
                                        }
                                    }
                                } else {
                                    PriorityQueue<ResultEntry> temp = new PriorityQueue<>(rankNum * 2);
                                    temp.add(new ResultEntry(score, peptide, true, true));
                                    ptmOnlyResultMap.put(scanNum, temp);
                                }
                            }
                        }
                    }
                }
            }
        }

        if (PIPI.DEV) {
            Set<Integer> tempSet = new HashSet<>(ptmFreeResultMap.keySet());
            tempSet.addAll(ptmOnlyResultMap.keySet());
            try (BufferedWriter writer = new BufferedWriter(new FileWriter("PTMFreeContainingGlobalScoreBound" + "." + batchStartIdx + "." + Thread.currentThread().getId() + ".csv"))) {
                writer.write("scan_num,PTM_free_low,PTM_free_high,PTM_containing_low,PTM_containing_high\n");
                for (int scanNum : tempSet) {
                    writer.write(scanNum + ",");
                    if (ptmFreeResultMap.containsKey(scanNum)) {
                        double highScore = 0;
                        double lowScore = 1;
                        for (ResultEntry temp : ptmFreeResultMap.get(scanNum)) {
                            if (temp.score > highScore) {
                                highScore = temp.score;
                            }
                            if (temp.score < lowScore) {
                                lowScore = temp.score;
                            }
                        }
                        writer.write(String.format("%f,%f,", lowScore, highScore));
                    } else {
                        writer.write("-,-,");
                    }
                    if (ptmOnlyResultMap.containsKey(scanNum)) {
                        double highScore = 0;
                        double lowScore = 1;
                        for (ResultEntry temp : ptmOnlyResultMap.get(scanNum)) {
                            if (temp.score > highScore) {
                                highScore = temp.score;
                            }
                            if (temp.score < lowScore) {
                                lowScore = temp.score;
                            }
                        }
                        writer.write(String.format("%f,%f\n", lowScore, highScore));
                    } else {
                        writer.write("-,-\n");
                    }
                }
            } catch (IOException ex) {
                ex.printStackTrace();
                System.exit(1);
            }
        }

        ptmFreeResult = convertResult(ptmFreeResultMap, peptideProteinMap, proteinSeqMap);
        ptmOnlyResult = convertResult(ptmOnlyResultMap, peptideProteinMap, proteinSeqMap);

        if (PIPI.DEV) {
            logger.debug("Writing searching results...");
            writeResult("ptm_only_candidates" + "." + batchStartIdx + "." + Thread.currentThread().getId() + ".csv", ptmOnlyResult);
        }
    }

    private void writeResult(String filePath,  Map<Integer, List<Peptide>> ptmOnly) {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(filePath))) {
            writer.write("scanNum,peptide,globalRank,is_decoy\n");
            for (int scanNum : ptmOnly.keySet()) {
                List<Peptide> peptideList = ptmOnly.get(scanNum);
                for (Peptide peptide : peptideList) {
                    if (peptide.isDecoy()) {
                        writer.write(scanNum + "," + peptide.getPTMFreeSeq() + "," + peptide.getGlobalRank() + ",1\n");
                    } else {
                        writer.write(scanNum + "," + peptide.getPTMFreeSeq() + "," + peptide.getGlobalRank() + ",0\n");
                    }
                }
            }
        } catch (IOException ex) {
            ex.printStackTrace();
            logger.error(ex.getMessage());
            System.exit(1);
        }
    }

    private Map<Integer, List<Peptide>> convertResult(Map<Integer, PriorityQueue<ResultEntry>> inputMap, Map<String, Set<String>> peptideProteinMap, Map<String, String> proteinSeqMap) {
        Map<Integer, List<Peptide>> outputMap = new HashMap<>();
        for (int scanNum : inputMap.keySet()) {
            List<Peptide> peptideList = new LinkedList<>();
            PriorityQueue<ResultEntry> tempQueue = inputMap.get(scanNum);
            int globalRank = tempQueue.size();
            while (!tempQueue.isEmpty()) {
                ResultEntry temp = tempQueue.poll();
                String leftFlank = "-";
                String rightFlank = "-";
                String peptideString = temp.peptide;
                if (peptideProteinMap.containsKey(peptideString)) {
                    String proSeq = proteinSeqMap.get(peptideProteinMap.get(peptideString).iterator().next());
                    int startIdx = proSeq.indexOf(peptideString);
                    if (startIdx == -1) {
                        logger.debug("Something wrong happened in Search.java (line: 223), scan num = {}, peptide = {}.", scanNum, peptideString);
                    } else if (startIdx == 0) {
                        int tempIdx = peptideString.length();
                        if (tempIdx < proSeq.length()) {
                            rightFlank = proSeq.substring(tempIdx, tempIdx + 1);
                        } else {
                            rightFlank = "-";
                        }
                    } else if (startIdx == proSeq.length() - peptideString.length()) {
                        leftFlank = proSeq.substring(startIdx - 1, startIdx);
                    } else {
                        leftFlank = proSeq.substring(startIdx - 1, startIdx);
                        rightFlank = proSeq.substring(startIdx + peptideString.length(), startIdx + peptideString.length() + 1);
                    }
                }
                peptideList.add(new Peptide(temp.peptide, temp.isDecoy(), massToolObj, maxMs2Charge, temp.score, leftFlank, rightFlank, globalRank, temp.containPTM()));
                --globalRank;
            }
            outputMap.put(scanNum, peptideList);
        }
        return outputMap;
    }

    public Map<Integer, List<Peptide>> getPTMOnlyResult() {
        return ptmOnlyResult;
    }

    public Map<Integer, List<Peptide>> getPTMFreeResult() {
        return ptmFreeResult;
    }
    public Map<Integer, Double> getTruthScore() {
        return truthScore;
    }

    public Set<Integer> getScanSearched() {return scanSearched;}
}
