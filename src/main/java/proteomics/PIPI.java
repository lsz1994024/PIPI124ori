package proteomics;

import ProteomicsLibrary.Types.PepWithScore;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.FDR.EstimateFDR;
import proteomics.Segment.InferenceSegment;
import proteomics.Spectrum.FilterSpectra;
import proteomics.Types.*;
import proteomics.Index.BuildIndex;
import proteomics.Parameter.Parameter;
import proteomics.Spectrum.PreSpectra;
import proteomics.TheoSeq.MassTool;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;
import uk.ac.ebi.pride.tools.mgf_parser.MgfFile;
import uk.ac.ebi.pride.tools.mzxml_parser.MzXMLFile;
import uk.ac.ebi.pride.tools.mzxml_parser.MzXMLParsingException;

import java.io.*;
import java.sql.*;
import java.text.DecimalFormat;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

public class PIPI {

    private static final Logger logger = LoggerFactory.getLogger(PIPI.class);
    private static final float C13Diff = 1.00335483f;
    private static final String versionStr = "1.2.4";

    public static final boolean DEV = false;

    public static void main(String args[]) {
        // Process inputs
        if (args.length != 2) {
            help();
        }

        // Set parameters
        String parameterPath = args[0].trim();
        String spectraPath = args[1].trim();

        try {
            logger.info("Running PIPI version {}.", versionStr);
            logger.info("Spectra: {}, parameter: {}", spectraPath, parameterPath);

            if (DEV) {
                logger.info("In DEV mode.");
            }

            new PIPI(parameterPath, spectraPath);
        } catch (Exception ex) {
            ex.printStackTrace();
            logger.error(ex.getMessage());
            System.exit(1);
        }
    }

    private PIPI(String parameterPath, String spectraPath) throws Exception {
        // Get the parameter map
        Parameter parameter = new Parameter(parameterPath);
        Map<String, String> parameterMap = parameter.returnParameterMap();
        float ms2Tolerance = Float.valueOf(parameterMap.get("ms2_tolerance"));
        float ms1Tolerance = Float.valueOf(parameterMap.get("ms1_tolerance"));
        int ms1ToleranceUnit = Integer.valueOf(parameterMap.get("ms1_tolerance_unit"));
        int maxMs2Charge = Integer.valueOf(parameterMap.get("max_ms2_charge"));
        int batchSize = Integer.valueOf(parameterMap.get("batch_size"));
        String percolatorPath = parameterMap.get("percolator_path");
        boolean outputPercolatorInput = (DEV || (Integer.valueOf(parameterMap.get("output_percolator_input")) == 1));

        logger.info("Indexing protein database...");
        BuildIndex buildIndexObj = new BuildIndex(parameterMap);
        MassTool massToolObj = buildIndexObj.returnMassToolObj();

        logger.info("Reading PTM database...");
        Map<String, TreeSet<Integer>> siteMass1000Map = readPTMDb(parameterMap, buildIndexObj.returnFixModMap());

        int tempMax = 0;
        int tempMin = 99999;
        for (TreeSet<Integer> tempSet : siteMass1000Map.values()) {
            if (tempSet.last() > tempMax) {
                tempMax = tempSet.last();
            }
            if (tempSet.first() < tempMin) {
                tempMin = tempSet.first();
            }
        }

        float minPtmMass = Float.valueOf(parameterMap.get("min_ptm_mass"));
        float maxPtmMass = Float.valueOf(parameterMap.get("max_ptm_mass"));

        logger.info("Reading spectra...");
        JMzReader spectraParser = null;
        String ext = "";
        try {
            File spectraFile = new File(spectraPath);
            if ((!spectraFile.exists() || (spectraFile.isDirectory()))) {
                throw new FileNotFoundException("The spectra file not found.");
            }
            String[] temp = spectraPath.split("\\.");
            ext = temp[temp.length - 1];
            if (ext.contentEquals("mzXML")) {
                spectraParser = new MzXMLFile(spectraFile);
            } else if (ext.contentEquals("mgf")) {
                spectraParser = new MgfFile(spectraFile);
            }
        } catch (FileNotFoundException | MzXMLParsingException ex) {
            ex.printStackTrace();
            logger.error(ex.getMessage());
            System.exit(1);
        }

        FilterSpectra filterSpectraObj = new FilterSpectra(spectraParser, parameterMap, massToolObj.returnMassTable());
        FilterSpectra.MassScan[] scanNumArray = filterSpectraObj.getScanNumArray();

        logger.info("Useful MS/MS spectra number: {}", scanNumArray.length);

        BufferedReader parameterReader = new BufferedReader(new FileReader("/home/slaiad/Data/PXD001468/9330/Truth9267.txt"));
//        BufferedReader parameterReader = new BufferedReader(new FileReader("/home/slaiad/Data/PXD022999/Truth.txt"));
        Map<Integer, String> pepTruth = new HashMap<>();
        Map<Integer, String> newTruth = new HashMap<>();
        String line;
        while ((line = parameterReader.readLine()) != null) {
            line = line.trim();
            String[] splitRes = line.split(",");
            pepTruth.put(Integer.valueOf(splitRes[0]), splitRes[1].replace('L','I'));
//            newTruth.put(Integer.valueOf(splitRes[0]), splitRes[1].replace('L','I'));
        }
        logger.info("Truth Loaded");

        int threadNum = Integer.valueOf(parameterMap.get("thread_num"));
        if (threadNum == 0) {
            threadNum = 1 + Runtime.getRuntime().availableProcessors();
        }
        ExecutorService threadPool = Executors.newFixedThreadPool(threadNum);

        InferenceSegment inference3SegmentObj = new InferenceSegment(buildIndexObj, ms2Tolerance, 3);

        List<FinalResultEntry> finalScoredPsms = new LinkedList<>();
        Map<Integer, ChargeMassTuple> numChargeMassMap = new HashMap<>();
        int startIdx;
        int endIdx = 0;
//        System.out.println(" LSZ here System");
        while (endIdx < scanNumArray.length) {
            startIdx = endIdx;
            if (batchSize == 0) {
                endIdx = scanNumArray.length;
            } else {
                endIdx = Math.min(startIdx + batchSize, scanNumArray.length);
            }

            logger.info("Searching batch {} - {} ({}%)...", startIdx, endIdx, String.format("%.1f", (float) endIdx * 100 / (float) scanNumArray.length));
//            for (FilterSpectra.MassScan scan : scanNumArray){
//            }
            PreSpectra preSpectraObj = new PreSpectra(spectraParser, parameterMap, massToolObj, ext, scanNumArray, startIdx, endIdx);
            Map<Integer, SpectrumEntry> numSpectrumMap = preSpectraObj.returnNumSpectrumMap();
            TreeMap<Float, List<Integer>> massNumMap = preSpectraObj.returnMassNumMap();
            numChargeMassMap.putAll(preSpectraObj.getNumChargeMassMap());

            if (!massNumMap.isEmpty()) {
                int batchSize2 = (massNumMap.size() / threadNum) + 1;
                Float[] massArray = massNumMap.keySet().toArray(new Float[massNumMap.size()]);
                Collection<PIPIWrap> taskList = new LinkedList<>();
                for (int i = 0; i < threadNum; ++i) {
                    int leftIdx = i * batchSize2;
                    int rightIdx = Math.min((i + 1) * batchSize2, massArray.length - 1);
                    if (leftIdx > massArray.length - 1) {
                        break;
                    }

                    if (rightIdx < massArray.length - 1) {
                        NavigableMap<Float, List<Integer>> subMassNumMap = massNumMap.subMap(massArray[leftIdx], true, massArray[rightIdx], false);
                        taskList.add(new PIPIWrap(buildIndexObj, massToolObj, inference3SegmentObj, numSpectrumMap, subMassNumMap, siteMass1000Map, ms1Tolerance, ms1ToleranceUnit, ms2Tolerance, minPtmMass, maxPtmMass, maxMs2Charge, startIdx, pepTruth, newTruth));
                    } else {
                        NavigableMap<Float, List<Integer>> subMassNumMap = massNumMap.subMap(massArray[leftIdx], true, massArray[rightIdx], true);
                        taskList.add(new PIPIWrap(buildIndexObj, massToolObj, inference3SegmentObj, numSpectrumMap, subMassNumMap, siteMass1000Map, ms1Tolerance, ms1ToleranceUnit, ms2Tolerance, minPtmMass, maxPtmMass, maxMs2Charge, startIdx, pepTruth, newTruth));
                    }
                }
                try {
                    List<Future<List<FinalResultEntry>>> tempResultList = threadPool.invokeAll(taskList);
                    for (Future<List<FinalResultEntry>> tempResult : tempResultList) {
                        if (tempResult.isDone() && !tempResult.isCancelled()) {
                            finalScoredPsms.addAll(tempResult.get());
                        } else {
                            logger.error("Threads were not finished normally.");
                            System.exit(1);
                        }
                    }
                } catch (Exception ex) {
                    ex.printStackTrace();
                    logger.error(ex.getMessage());
                    System.exit(1);
                }
            }
        }

        // shutdown threads.
        threadPool.shutdown();
        try {
            if (!threadPool.awaitTermination(60, TimeUnit.SECONDS)) {
                threadPool.shutdownNow();
                if (!threadPool.awaitTermination(60, TimeUnit.SECONDS))
                    System.err.println("Pool did not terminate");
            }
        } catch (InterruptedException ie) {
            threadPool.shutdownNow();
            Thread.currentThread().interrupt();
            logger.error("Threads were not finished normally.");
            System.exit(1);
        }

        logger.info("Estimating FDR...");
        // estimate T-D FDR
        PFM(finalScoredPsms, buildIndexObj.returnPepProMap(), buildIndexObj.protPepNum, pepTruth, "9268_TD");
        new EstimateFDR(finalScoredPsms);

        // estimate Percolator FDR
        String percolatorInputFileName = spectraPath + ".input.temp";
        String percolatorOutputFileName = spectraPath + ".output.temp";
        Map<String, Set<String>> peptideProteinMap = buildIndexObj.returnPepProMap();
        Map<String, String> decoyPeptideProteinMap = buildIndexObj.returnDecoyPepProMap();
//        writePercolator(finalScoredPsms, numChargeMassMap, peptideProteinMap, decoyPeptideProteinMap, percolatorInputFileName, buildIndexObj.returnFixModMap());
        Map<Integer, PercolatorEntry> percolatorResultMap = runPercolator(percolatorPath, percolatorInputFileName, percolatorOutputFileName);

        if (percolatorResultMap.isEmpty()) {
            logger.info("Percolator failed to estimate FDR. The results won't contain percolator_score, posterior_error_prob, and percolator_q_value.");
        }

        if (!outputPercolatorInput) {
            (new File(percolatorInputFileName)).delete();
            (new File(percolatorOutputFileName)).delete();
        }

        logger.info("Saving results...");
        writeFinalResult(finalScoredPsms, percolatorResultMap, numChargeMassMap, peptideProteinMap, spectraPath + ".pipi.csv", buildIndexObj.returnFixModMap());

        logger.info("Done.");
    }

    private void PFM(List<FinalResultEntry> resultList, Map<String, Set<String>> peptide0Map,  Map<String, Integer> protPepNum, Map<Integer, String> pepTruth, String srchId) {

        System.out.println("PFM starts ========================");
        Map<String, Integer> pepTruthTimes = new HashMap<>();
        for (int scanNum : pepTruth.keySet()){
            String pepSeq = pepTruth.get(scanNum).replace('L', 'I');
            if (pepTruthTimes.containsKey(pepSeq)){
                pepTruthTimes.put(pepSeq, pepTruthTimes.get(pepSeq)+1);
            }else {
                pepTruthTimes.put(pepSeq, 1);
            }
        }
        Multimap<String,String> protIdInTruth = HashMultimap.create();

        for(String pep : pepTruthTimes.keySet()){
//            int pepTimes = pepTruthTimes.get(pep);
//            System.out.println("pepTimes " + pep +" " + pepTimes );
            if (!peptide0Map.containsKey(pep)){
                System.out.println("not in 1 "+pep );
                continue;
            }
            for (String proId : peptide0Map.get(pep)) {
                protIdInTruth.put(proId,pep);
            }
        }

//        for (String prot : protIdInTruth.keySet()) {
//            System.out.println(prot.substring(0,10)+": "+protIdInTruth.get(prot).size() + ", "+protPepNum.get(prot));
//            for (String pep : protIdInTruth.get(prot)){
//                System.out.println(pep + ", pepTime: " + pepTruthTimes.get(pep) + ", pepProt Num: " + pepProtMap.get(pep).size());
//            }
//        }
        System.out.println("====== truth up ========= real down ============");
        System.out.println("num of result psms "+ resultList.size());
        TreeMap<Integer, String> candsStrs = new TreeMap<>();
        List<PepWithScore> pepWithScoreList = new LinkedList<>();
        for (FinalResultEntry psm : resultList){
            int scanNum = psm.getScanNum();
            List<PepWithScore> pepCandisList = psm.pepCandisList;
            Collections.sort(pepCandisList);

            String tempStr = scanNum+",";
            for (PepWithScore pepcand : pepCandisList){
//                System.out.println(scanNum);
                tempStr += pepcand.pepSeq+","+pepcand.score+","+pepcand.isDecoy+","+ pepcand.hasPTM+","+String.join("-",pepcand.proteins)+",";
            }
            candsStrs.put(scanNum, tempStr.substring(0, tempStr.length()-1)+"\n");
            List<PepWithScore> tempBackbonesList = new LinkedList<>();
            tempBackbonesList.add(pepCandisList.get(0));
            double tempTopScore = pepCandisList.get(0).score;
            for (PepWithScore pep : pepCandisList.subList(1, pepCandisList.size()-1)) {
                if (pep.score < tempTopScore){ break;}
                tempBackbonesList.add(pep);
            }
            double count = 1.0/tempBackbonesList.size();
            for (PepWithScore pep : tempBackbonesList){
//                pepWithScoreList.add(new PepWithScore(pep.pepSeq,tempTopScore,count, false));

            }
        }

        System.out.println("Start writing ========================");
        String resultPath = srchId + ".candis.csv";
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(resultPath))) {
            writer.write("scanNo,pep1,s1,b1,m1,p1,pep2,s2,b2,m2,p2,pep3,s3,b3,m3,p3,pep4,s4,b4,m4,p4,pep5,s5,b5,m5,p5,pep6,s6,b6,m6,p6,pep7,s7,b7,m7,p7,pep8,s8,b8,m8,p8,pep9,s9,b9,m9,p9,pep10,s10,b10,m10,p10,pep11,s11,b11,m11,p11,pep12,s12,b12,m12,p12,pep13,s13,b13,m13,p13,pep14,s14,b14,m14,p14,pep15,s15,b15,m15,p15,pep16,s16,b16,m16,p16,pep17,s17,b17,m17,p17,pep18,s18,b18,m18,p18,pep19,s19,b19,m19,p19,pep20,s20,b20,m20,p20\n");
            for (int scanNo : candsStrs.keySet()) {
                writer.write(candsStrs.get(scanNo));
            }
        } catch (IOException ex) {
            ex.printStackTrace();
            logger.error(ex.getMessage());
            System.exit(1);
        }
        System.out.println("End writing ========================");
        Collections.sort(pepWithScoreList);

        Map<String, Double> protScoreNewMap = new HashMap<>();
        for (PepWithScore pep : pepWithScoreList.subList(0, (int) Math.ceil(pepWithScoreList.size() * 0.1))) {
            String pepSeq = pep.pepSeq;
            int numPepProt = peptide0Map.get(pepSeq).size();
            double countScore = pep.score/numPepProt*1;//*pep.count
            DecimalFormat df = new DecimalFormat("#.####");
            double roundedScore = Double.valueOf(df.format(countScore));
            for (String proId : peptide0Map.get(pepSeq)) {
//                System.out.println(pepSeq+" count: " + pepCountMap.get(pepSeq)+ " prot: "+proId+" , ");
                if (protScoreNewMap.containsKey(proId)) {
                    protScoreNewMap.put(proId, protScoreNewMap.get(proId) + roundedScore);
                } else {
                    protScoreNewMap.put(proId, roundedScore);
                }
            }
        }
        System.out.println("the last top-hit considered "+ pepWithScoreList.get((int) Math.ceil(pepWithScoreList.size() * 0.1)).score);

        for (String prot : protScoreNewMap.keySet()){
            System.out.println(prot+" : "+ protScoreNewMap.get(prot));
        }


        Map<String, Double> pepCountMap = new HashMap<>();
        for (PepWithScore pep : pepWithScoreList.subList(0, (int) Math.ceil(pepWithScoreList.size() * 0.02))) {
            if (pepCountMap.containsKey(pep.pepSeq)){
                pepCountMap.put(pep.pepSeq, pepCountMap.get(pep.pepSeq)+pep.count);
            } else {
                pepCountMap.put(pep.pepSeq, pep.count);
            }
//            System.out.println(pep.pepSeq + " : "+pep.score+ " , "+pep.count);
        }

        int maxPepNumOneProt = 0;
        String maxProtId = "";

        for (String pep : pepCountMap.keySet()) {
//            System.out.println("pepCount" + pep+" , "+ pepCountMap.get(pep));
            for (String proId : peptide0Map.get(pep)) {
                proId = proId.trim();
                if (protPepNum.get(proId) > maxPepNumOneProt){
                    maxPepNumOneProt = protPepNum.get(proId);
                    maxProtId = proId;
                }
            }
        }
        System.out.println("maxPepNumOneProt "+maxPepNumOneProt);
        maxPepNumOneProt += 100;

        Map<String, Double> protScoreMap = new HashMap<>();
        Multimap<String,String> protPepMap = HashMultimap.create();
        for (String pepSeq : pepCountMap.keySet()){
            for (String proId : peptide0Map.get(pepSeq)) {
                double countScore = countScore(pepCountMap.get(pepSeq),protPepNum.get(proId), maxPepNumOneProt);
                DecimalFormat df = new DecimalFormat("#.####");
                double roundedScore = Double.valueOf(df.format(countScore));
                System.out.println(pepSeq+" count: " + pepCountMap.get(pepSeq)+ " prot: "+proId+" , ");
                if (protScoreMap.containsKey(proId)) {
                    protScoreMap.put(proId, protScoreMap.get(proId) + roundedScore);
                } else {
                    protScoreMap.put(proId, roundedScore);
                }
                protPepMap.put(proId, pepSeq);
            }
        }


//        int truthProtCovered = 0;
//        int wrongProtImported = 0;
        for (String prot : protScoreMap.keySet()) {
//            int isTrueProt;
//            if (protIdInTruth.containsKey(prot)){
//                truthProtCovered ++;
//                isTrueProt = 1;
//            }else{
//                wrongProtImported ++;
//                isTrueProt = 0;
//            }
            System.out.println(prot + ": "+protScoreMap.get(prot));
            for (String pep : protPepMap.get(prot)){
                System.out.println(pep + ": "+pepCountMap.get(pep));
            }
        }
//        System.out.println("truthProtCovered "+truthProtCovered+ " wrongProtImported "+ wrongProtImported);
        System.out.println("protScoreMap "+protScoreMap.size());

        for (FinalResultEntry psm : resultList){
            int scanNum = psm.getScanNum();

            if (!pepTruth.containsKey(scanNum)){ // only calculate accuracy for scans in simulation 1
//                System.out.println("not in 2" );
                continue;
            }
            if (scanNum == 13701 || scanNum == 2095 || scanNum == 2052 || scanNum == 22598){
                System.out.println("here");
            }
            List<PepWithScore> pepCandisList = psm.pepCandisList;

            List<PepWithScore> pepWithProtScore = new LinkedList<>();

            for (PepWithScore pep : pepCandisList) {
                if (!peptide0Map.containsKey(pep.pepSeq)){
                    System.out.println("not in 3");
                    continue;
                }
                String pepSeq = pep.pepSeq;
                if (pepSeq.equals(pepTruth.get(scanNum))){
                    System.out.println(scanNum+" truthIn20");
                }
                double bestProtScore = -1d;
                for (String prot : peptide0Map.get(pepSeq)) {
                    prot = prot.trim();
                    double protScore = -1d;
                    if (protScoreNewMap.containsKey(prot)){
                        protScore = protScoreNewMap.get(prot);
                    }
                    if (protScore > bestProtScore){
                        bestProtScore = protScore;
                    }
                }
//                pepWithProtScore.add(new PepWithScore(pepSeq, bestProtScore, pep.score, false));
            }

            Collections.sort(pepWithProtScore);
            double bestScore = pepWithProtScore.get(0).score;
            String bestSeq = pepWithProtScore.get(0).pepSeq;
            double bestPepScore = pepWithProtScore.get(0).count;
            for (PepWithScore pep : pepWithProtScore){
                String pepSeq = pep.pepSeq;
                double score = pep.score;
                double pepScore = pep.count;
                if  (score == bestScore) {
                    if (pepScore > bestPepScore) {
                        bestSeq = pepSeq;
                        bestPepScore = pepScore;
                    }
                } else {
                    break;
                }
            }

//            System.out.println(" Here");
            if (bestSeq.isEmpty()){
                System.out.println(" empty continue" +  scanNum);
                continue;
//                bestSeq = oriBestBackbone;
//                if (bestSeq.isEmpty()){continue;}
            }
            if (bestSeq.equals(pepTruth.get(scanNum))){
                System.out.println(scanNum+" yes");
            }else{
                System.out.println(scanNum+" no");
            }
        }
        System.out.println("Finished PFM================");
    }

    private double countScore(double count, int k, int bound){
        if (bound<=k){
            System.out.println("bound<=k " + bound+ " "+ k);
        }
        double ei = Math.exp(count);
        double elog = Math.exp( (Math.log(((double)bound-k)/(bound+k)) + 1) * count );
        return bound*(ei-elog)/(k*(ei+elog));
    }

    private static void help() {
        String helpStr = "PIPI version " + versionStr + "\r\n"
                + "A tool identifying peptides with unlimited PTM.\r\n"
                + "Author: Fengchao Yu\r\n"
                + "Email: fyuab@connect.ust.hk\r\n"
                + "ECL usage: java -Xmx25g -jar /path/to/PIPI.jar <parameter_file> <data_file>\r\n"
                + "\t<parameter_file>: parameter file. Can be download along with PIPI.\r\n"
                + "\t<data_file>: spectra data file (mzXML)\r\n"
                + "\texample: java -Xmx32g -jar ECL.jar parameter.def data.mzxml\r\n";
        System.out.print(helpStr);
        System.exit(1);
    }

    private static void writePercolator(List<FinalResultEntry> finalScoredResult, Map<Integer, ChargeMassTuple> numChargeMassMap, Map<String, Set<String>> peptideProteinMap, Map<String, String> decoyPeptideProteinMap, String resultPath, Map<String, Float> fixModMap) {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(resultPath))) {
            writer.write("id\tlabel\tscannr\txcorr\tdelta_c\tdelta_L_c\tnegative_log10_e_value\tnormalized_cross_corr\tglobal_search_rank\tabs_ppm\tIonFrac\tmatched_high_peak_frac\tcharge1\tcharge2\tcharge3\tcharge4\tcharge5\tcharge6\tvar_PTM_num\tpeptide\tprotein\n");
            for (FinalResultEntry entry : finalScoredResult) {
                float expMass = numChargeMassMap.get(entry.getScanNum()).mass;
                Peptide peptide = entry.getPeptide();
                float theoMass = peptide.getPrecursorMass();
                float massDiff = getMassDiff(expMass, theoMass, C13Diff);
                String proteinIdStr = "";

                if (!entry.isDecoy()) {
                    for (String temp : peptideProteinMap.get(peptide.getPTMFreeSeq())) {
                        proteinIdStr += temp + ";";
                    }
                } else {
                    proteinIdStr = decoyPeptideProteinMap.get(peptide.getPTMFreeSeq());
                }

                StringBuilder sb = new StringBuilder(20);
                int charge = numChargeMassMap.get(entry.getScanNum()).charge;
                for (int i = 0; i < 6; ++i) {
                    if (i == charge - 1) {
                        sb.append(1);
                    } else {
                        sb.append(0);
                    }
                    sb.append("\t");
                }

                if (entry.isDecoy()) {
                    writer.write(entry.getScanNum() + "\t-1\t" + entry.getScanNum() + "\t" + entry.getScore() + "\t" + entry.getDeltaC() + "\t" + entry.getDeltaLC() + "\t" + entry.getNegativeLog10EValue() + "\t" + entry.getNormalizedCrossXcorr() + "\t" + entry.getGlobalSearchRank() + "\t" + Math.abs(massDiff * 1e6f / theoMass) + "\t" + entry.getIonFrac() + "\t" + entry.getMatchedHighestIntensityFrac() + "\t" + sb.toString() + peptide.getVarPTMNum() + "\t" + peptide.getLeftFlank() + "." + peptide.getPTMContainedString(fixModMap) + "." + peptide.getRightFlank() + "\t" + proteinIdStr + "\n");
                } else {
                    writer.write(entry.getScanNum() + "\t1\t" + entry.getScanNum() + "\t" + entry.getScore() + "\t" + entry.getDeltaC() + "\t" + entry.getDeltaLC() + "\t" + entry.getNegativeLog10EValue() + "\t" + entry.getNormalizedCrossXcorr() + "\t" + entry.getGlobalSearchRank() + "\t" + Math.abs(massDiff * 1e6f / theoMass) + "\t" + entry.getIonFrac() + "\t" + entry.getMatchedHighestIntensityFrac() + "\t" + sb.toString() + peptide.getVarPTMNum() + "\t" + peptide.getLeftFlank() + "." + peptide.getPTMContainedString(fixModMap) + "." + peptide.getRightFlank() + "\t" + proteinIdStr + "\n");
                }
            }
        } catch (IOException ex) {
            ex.printStackTrace();
            logger.error(ex.getMessage());
            System.exit(1);
        }
    }

    private static Map<Integer, PercolatorEntry> runPercolator(String percolatorPath, String percolatorInputFileName, String percolatorOutputFileName) {
        Map<Integer, PercolatorEntry> percolatorResultMap = new HashMap<>();
        try {
            if ((new File(percolatorPath)).exists()) {
                Process ps = Runtime.getRuntime().exec(percolatorPath + " --only-psms --results-psms " + percolatorOutputFileName + " " + percolatorInputFileName);
                ps.waitFor();

                if (!(new File(percolatorOutputFileName).exists())) {
                    logger.warn("Error in running Percolator. The results won't contain percolator_score, posterior_error_prob, and percolator_q_value.");
                    return percolatorResultMap;
                }

                BufferedReader reader = new BufferedReader(new FileReader(percolatorOutputFileName));
                String line;
                while ((line = reader.readLine()) != null) {
                    line = line.trim();
                    if (!line.startsWith("PSMId")) {
                        String[] parts = line.split("\t");
                        percolatorResultMap.put(Integer.valueOf(parts[0]), new PercolatorEntry(Double.valueOf(parts[1]), parts[2], parts[3]));
                    }
                }
                reader.close();

                if (DEV) {
                    reader = new BufferedReader(new InputStreamReader(ps.getErrorStream()));
                    while ((line = reader.readLine()) != null) {
                        System.err.print(line);
                    }
                    reader.close();
                }
            } else {
                logger.error("Cannot find Percolator for estimating Percolator Q-Value. The results won't contain percolator_score, posterior_error_prob, and percolator_q_value.");
                return percolatorResultMap;
            }
        } catch (IOException | InterruptedException ex) {
            ex.printStackTrace();
            logger.error(ex.getMessage());
            return percolatorResultMap;
        }

        return percolatorResultMap;
    }

    private static void writeFinalResult(List<FinalResultEntry> finalScoredPsms, Map<Integer, PercolatorEntry> percolatorResultMap, Map<Integer, ChargeMassTuple> numChainMassMap, Map<String, Set<String>> peptideProteinMap, String outputPath, Map<String, Float> fixModMap) {
        TreeMap<Double, List<String>> tempMap = new TreeMap<>();
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputPath))) {
            writer.write("scan_num,peptide,charge,theo_mass,exp_mass,ppm,protein_ID,xcorr,e_value,naive_q_value,percolator_score,posterior_error_prob,percolator_q_value\n");
            for (FinalResultEntry entry : finalScoredPsms) {
                if (!entry.isDecoy()) {
                    int scanNum = entry.getScanNum();
                    float expMass = numChainMassMap.get(scanNum).mass;
                    int charge = numChainMassMap.get(scanNum).charge;
                    Peptide peptide = entry.getPeptide();
                    float theoMass = peptide.getPrecursorMass();
                    float massDiff = getMassDiff(expMass, theoMass, C13Diff);
                    float ppm = Math.abs(massDiff * 1e6f / theoMass);


                    String proteinIdStr = "";
                    for (String tempStr : peptideProteinMap.get(peptide.getPTMFreeSeq())) {
                        proteinIdStr += tempStr + ";";
                    }

                    String str;
                    boolean sortedByPercolatorScore = true;
                    if (percolatorResultMap.containsKey(scanNum)) {
                        PercolatorEntry percolatorEntry = percolatorResultMap.get(scanNum);
                        str = String.format("%d,%s,%d,%.4f,%.4f,%.2f,%s,%.4f,%E,%.4f,%.3f,%s,%s\n", scanNum, peptide.getPTMContainedString(fixModMap), charge, theoMass, expMass, ppm, proteinIdStr, entry.getScore(), entry.getEValue(), entry.getQValue(), percolatorEntry.percolatorScore, percolatorEntry.PEP, percolatorEntry.qValue);
                    } else {
                        sortedByPercolatorScore = false;
                        str = String.format("%d,%s,%d,%.4f,%.4f,%.2f,%s,%.4f,%E,%.4f,%s,%s,%s\n", scanNum, peptide.getPTMContainedString(fixModMap), charge, theoMass, expMass, ppm, proteinIdStr, entry.getScore(), entry.getEValue(), entry.getQValue() , "-", "-", "-");
                    }

                    if (sortedByPercolatorScore) {
                        if (tempMap.containsKey(percolatorResultMap.get(scanNum).percolatorScore)) {
                            tempMap.get(percolatorResultMap.get(scanNum).percolatorScore).add(str);
                        } else {
                            List<String> tempList = new LinkedList<>();
                            tempList.add(str);
                            tempMap.put(percolatorResultMap.get(scanNum).percolatorScore, tempList);
                        }
                    } else {
                        if (tempMap.containsKey(entry.getNegativeLog10EValue())) {
                            tempMap.get(entry.getNegativeLog10EValue()).add(str);
                        } else {
                            List<String> tempList = new LinkedList<>();
                            tempList.add(str);
                            tempMap.put(entry.getNegativeLog10EValue(), tempList);
                        }
                    }
                }
            }

            Double[] tempArray = tempMap.keySet().toArray(new Double[tempMap.size()]);
            for (int i = tempArray.length - 1; i >= 0; --i) {
                List<String> tempList = tempMap.get(tempArray[i]);
                for (String tempStr : tempList) {
                    writer.write(tempStr);
                }
            }
        } catch (IOException ex) {
            ex.printStackTrace();
            logger.error(ex.getMessage());
            System.exit(1);
        }
    }

    private static float getMassDiff(float expMass, float theoMass, float C13Diff) {
        float massDiff1 = expMass - theoMass;
        float massDiff2 = expMass - theoMass - C13Diff;
        float massDiff3 = expMass - theoMass - 2 * C13Diff;
        float absMassDiff1 = Math.abs(massDiff1);
        float absMassDiff2 = Math.abs(massDiff2);
        float absMassDiff3 = Math.abs(massDiff3);

        if ((absMassDiff1 <= absMassDiff2) && (absMassDiff1 <= absMassDiff2)) {
            return massDiff1;
        } else if ((absMassDiff2 <= absMassDiff1) && (absMassDiff2 <= absMassDiff3)) {
            return massDiff2;
        } else {
            return massDiff3;
        }
    }

    private Map<String, TreeSet<Integer>> readPTMDb(Map<String, String> parameterMap, Map<String, Float> fixModMap) {
        Map<String, TreeSet<Integer>> siteMass1000Map = new HashMap<>();
        float minPtmMass = Float.valueOf(parameterMap.get("min_ptm_mass"));
        float maxPtmMass = Float.valueOf(parameterMap.get("max_ptm_mass"));

        try (BufferedReader reader = new BufferedReader(new FileReader(parameterMap.get("PTM_db")))) {
            String line;
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if (line.startsWith("site") || line.startsWith("#")) {
                    continue;
                }

                String[] parts = line.split("\t");
                String site = parts[0];
                String position = parts[3];
                float mass = Float.valueOf(parts[2]);

                if ((mass > maxPtmMass) || mass < minPtmMass) {
                    continue;
                }

                int mass1000 = (int) Math.floor(mass * 1000);

                String siteString;
                if (site.contentEquals("N-term") && (Math.abs(fixModMap.get("n")) < 1e-6)) {
                    if (position.contentEquals("PROTEIN_N")) {
                        siteString = "PROTEIN_N";
                    } else {
                        siteString = "PEPTIDE_N";
                    }
                } else if (site.contentEquals("C-term")) {
                    if (position.contentEquals("PROTEIN_C")) {
                        siteString = "PROTEIN_C";
                    } else {
                        siteString = "PEPTIDE_C";
                    }
                } else if (Math.abs(fixModMap.get(site)) < 1e-6) { // fix modified amino acid cannot be modified again.
                    if (position.contentEquals("PROTEIN_N")) {
                        siteString = site + "-PROTEIN_N";
                    } else if (position.contentEquals("PROTEIN_C")) {
                        siteString = site + "-PROTEIN_C";
                    } else if (position.contentEquals("PEPTIDE_N")) {
                        siteString = site + "-PEPTIDE_N";
                    } else if (position.contentEquals("PEPTIDE_C")) {
                        siteString = site + "-PEPTIDE_C";
                    } else {
                        siteString = site;
                    }
                } else {
                    continue;
                }

                if (siteMass1000Map.containsKey(siteString)) {
                    siteMass1000Map.get(siteString).add(mass1000);
                } else {
                    TreeSet<Integer> tempSet = new TreeSet<>();
                    tempSet.add(mass1000);
                    siteMass1000Map.put(siteString, tempSet);
                }
            }
        } catch (IOException ex) {
            ex.printStackTrace();
            logger.error(ex.getMessage());
            System.exit(1);
        }

        return siteMass1000Map;
    }
}
