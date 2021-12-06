package ProteomicsLibrary;

//import ProteomicsLibrary.Types.SparseVector;
import ProteomicsLibrary.Types.SparseBooleanVector;

import java.util.*;

public class PrepareSpectrum {

    private static final double defaultIntensity = 1;
    private static final int xcorrOffset = 75;
    private static final double removePrecursorPeakTolerance = 1.5;

    private final MassTool massTool;

    public PrepareSpectrum(MassTool massTool) {
        this.massTool = massTool;
    }

    public TreeMap<Double, Double> preSpectrumTopNStyle (Map<Double, Double> inputPL, double precursorMass, int precursorCharge, double minClear, double maxClear, int topN) {
        TreeMap<Double, Double> outputPL = removeCertainPeaks(inputPL, precursorMass, precursorCharge, minClear, maxClear);
        if (outputPL.subMap(0d, precursorMass).isEmpty()) {
            return new TreeMap<>();
        } else {
//            System.out.println("raw");
//            for (Map.Entry<Double, Double> pl : inputPL.entrySet() )
//            {
//                System.out.println(pl.getKey()+","+pl.getValue());
//            }

//            System.out.println("====================");
            return topNStyleNormalization(sqrtPL(new TreeMap<>(outputPL.subMap(0d, precursorMass))), topN);  // keep the top N peaks in the window size of 100 m/z
        }
    }

    public TreeMap<Double, Double> NoiseCalculate (Map<Double, Double> inputPL, double precursorMass, int precursorCharge, double minClear, double maxClear, int topN) {
        TreeMap<Double, Double> outputPL = removeCertainPeaks(inputPL, precursorMass, precursorCharge, minClear, maxClear);
        if (outputPL.subMap(0d, precursorMass).isEmpty()) {
            return new TreeMap<>();
        } else {
            return NoiseMap(sqrtPL(new TreeMap<>(outputPL.subMap(0d, precursorMass))), topN);  // keep the top N peaks in the window size of 100 m/z
        }
    }

    public TreeMap<Double, Double> preSpectrumTopNStyle (Map<Double, Double> inputPL, double precursorMass, int topN) {
        TreeMap<Double, Double> outputPL = new TreeMap<>(inputPL);
        if (outputPL.subMap(0d, precursorMass).isEmpty()) {
            return new TreeMap<>();
        } else {
            return topNStyleNormalization(sqrtPL(new TreeMap<>(outputPL.subMap(0d, precursorMass))), topN);
        }
    }

    public SparseBooleanVector prepareXCorr(TreeMap<Double, Double> plMap, boolean flankingPeaks) {
        if (plMap.isEmpty()) {
            return new SparseBooleanVector();
        } else {

            double[] plArray = new double[massTool.mzToBin(plMap.lastKey()) + 1];
            for (double mz : plMap.keySet()) {
                if (Math.abs(plMap.get(mz)) > 1e-6) {
                    int idx = massTool.mzToBin(mz);
                    plArray[idx] = Math.max(plMap.get(mz), plArray[idx]);
                }
            }
            return prepareXCorr(plArray, flankingPeaks);
        }
    }

    public SparseBooleanVector prepareXCorr(double[] plArray, boolean flankingPeaks) {
        SparseBooleanVector xcorrPL = new SparseBooleanVector();
        int offsetRange = 2 * xcorrOffset + 1;
        double factor = 1 / (double) (offsetRange - 1);
        double mySum = 0;
        for (int i = 0; i < xcorrOffset; ++i) {
            mySum += plArray[i];
        }

        double[] tempArray = new double[plArray.length];
        for (int i = xcorrOffset; i < plArray.length + xcorrOffset; ++i) {
            if (i < plArray.length) {
                mySum += plArray[i];
            }
            if (i >= offsetRange) {
                mySum -= plArray[i - offsetRange];
            }
            tempArray[i - xcorrOffset] = (mySum - plArray[i - xcorrOffset]) * factor;
        }

        for (int i = 1; i < plArray.length; ++i) {
            double temp = plArray[i] - tempArray[i];
            if (flankingPeaks) {
                temp += (plArray[i - 1] - tempArray[i - 1]) * 0.5;
                if (i + 1 < plArray.length) {
                    temp += (plArray[i + 1] - tempArray[i + 1]) * 0.5;
                }
            }
            if (Math.abs(temp) > 1e-6) {
                xcorrPL.put(i, temp);
            }
        }

        return xcorrPL;
    }

    public SparseBooleanVector digitizePL(TreeMap<Double, Double> plMap) {
        SparseBooleanVector digitizedPL = new SparseBooleanVector();
        for (double mz : plMap.keySet()) {
            int idx = massTool.mzToBin(mz);
            if (Math.abs(plMap.get(mz)) > 1e-6) {
                digitizedPL.put(idx, Math.max(plMap.get(mz), digitizedPL.get(idx)));
            }
        }
        return digitizedPL;
    }


    // keep the top N peaks in the window size of 100 m/z
    public static TreeMap<Double, Double> topNStyleNormalization(TreeMap<Double, Double> inputPL, int localTopN) {
        if (inputPL.isEmpty()) {
            return new TreeMap<>();
        } else {

            TreeMap<Double, Double> outputPL = new TreeMap<>();
            double minMz = inputPL.firstKey();
            double maxMz = inputPL.lastKey();
            double leftMz = minMz;
            int sectionId = 0;
            while (leftMz < maxMz) {

                double rightMz = Math.min(leftMz + 100, maxMz);  // not divide the mz range into 10 section (as said in paper). actually, every 100Hz
                NavigableMap<Double, Double> subMap;
                if (rightMz < maxMz) {
                    subMap = inputPL.subMap(leftMz, true, rightMz, false);
                } else {
                    subMap = inputPL.subMap(leftMz, true, rightMz, true);
                }
                if (!subMap.isEmpty()) {
                    Double[] intensityArray = subMap.values().toArray(new Double[0]);
                    Arrays.sort(intensityArray, Comparator.reverseOrder());
                    double temp1 = defaultIntensity / intensityArray[0];  // use this as a factor, to normalize all
                    double temp2 = subMap.size() > localTopN ? intensityArray[localTopN] : 0;  // if subMap has less than topN peaks, let temp2 = 0, so that every peak will be bigger than temp2, thus take them all

//                    System.out.println(sectionId+" maxIntens  "+ temp1);
                    for (double mz : subMap.keySet()) {
                        if (subMap.get(mz) > temp2) {
                            outputPL.put(mz, subMap.get(mz) * temp1);
//                            System.out.println(sectionId+" add "+mz+ " inte "+subMap.get(mz) * temp1);
                        }
                    }
                }
                leftMz = rightMz;
                sectionId ++;
            }

            return outputPL;
        }
    }

    public static TreeMap<Double, Double> NoiseMap(TreeMap<Double, Double> inputPL, int localTopN) {
        if (inputPL.isEmpty()) {
            return new TreeMap<>();
        }
        else {
            TreeMap<Double, Double> outputPL = new TreeMap<>();
            double minMz = inputPL.firstKey();
            double maxMz = inputPL.lastKey();
            double leftMz = minMz;
            while (leftMz < maxMz) {
                double rightMz = Math.min(leftMz + 100, maxMz);
                NavigableMap<Double, Double> subMap;
                if (rightMz < maxMz) {
                    subMap = inputPL.subMap(leftMz, true, rightMz, false);
                } else {
                    subMap = inputPL.subMap(leftMz, true, rightMz, true);
                }
                if (!subMap.isEmpty()) {
                    Double[] intensityArray = subMap.values().toArray(new Double[0]);
                    Arrays.sort(intensityArray, Comparator.reverseOrder());
                    double temp1 = defaultIntensity / intensityArray[0];
                    double temp2 = 0.0;
                    if (subMap.size() > localTopN + 1) {
                        temp2 = intensityArray[localTopN];
                        for (double mz : subMap.keySet()) {
                            if (subMap.get(mz) < temp2) {
                                Double tempSNR = 0.0;
                                if (outputPL.containsKey(leftMz)){
                                    tempSNR = subMap.get(mz) * temp1 + outputPL.get(leftMz);
                                    outputPL.put(leftMz, tempSNR);
                                }
                                else {
                                    tempSNR = subMap.get(mz) * temp1;
                                    outputPL.put(leftMz, tempSNR);
                                }
                            }
                        }
                        outputPL.put(leftMz, outputPL.get(leftMz) / (subMap.size() - localTopN));
//                        outputPL.put(leftMz, outputPL.get(leftMz));
                    }
                    else{
                        outputPL.put(leftMz, 1.0);
                    }
                }
                leftMz = rightMz;
            }

            return outputPL;
        }
    }

    private static TreeMap<Double, Double> removeCertainPeaks(Map<Double, Double> peakMap, double precursorMass, int precursorCharge, double minClear, double maxClear) {
        TreeMap<Double, Double> mzIntensityMap = new TreeMap<>();
        double precursorMz = precursorMass / precursorCharge + MassTool.PROTON;
//        System.out.println("precursorMass "+precursorMass);
//        System.out.println("precursorMz "+precursorMz);
        for (double mz : peakMap.keySet()) {
            if (((mz < minClear) || (mz > maxClear)) && (mz > 50)) {
                if ((peakMap.get(mz) > 1e-6) && (Math.abs(peakMap.get(mz) - precursorMz) > removePrecursorPeakTolerance)) {
                    mzIntensityMap.put(mz, peakMap.get(mz));  // when the peak is qualified, it is remained, otherwise, removed
                }
            }
        }

        return mzIntensityMap;
    }

    private static TreeMap<Double, Double> sqrtPL(TreeMap<Double, Double> plMap) {

        TreeMap<Double, Double> sqrtPlMap = new TreeMap<>();
        for (double mz : plMap.keySet()) {
            if (plMap.get(mz) > 1e-6) {
                sqrtPlMap.put(mz, Math.sqrt(plMap.get(mz)));
            }
        }
        return sqrtPlMap;
    }
}
