package ProteomicsLibrary;

import ProteomicsLibrary.Types.Coordinate;

import java.util.*;

public class Score {

    public static double calIonFraction(double[][] ionMatrix, int precursorCharge, Map<Double, Double> plMap, double ms2Tolerance) {
        int matchedPeakNum = 0;
        int maxRow = Math.max(2, Math.min(ionMatrix.length, 2 * (precursorCharge - 1)));
        int totalIonNum = ionMatrix[0].length * maxRow;
        for (int i = 0; i < maxRow; ++i) {
            for (int j = 0; j < ionMatrix[0].length; ++j) {
                for (double mz : plMap.keySet()) {
                    if (Math.abs(mz - ionMatrix[i][j]) <= ms2Tolerance) {
                        ++matchedPeakNum;
                        break;
                    }
                }
            }
        }

        return (double) matchedPeakNum / (double) totalIonNum;
    }

    public static double calMatchedHighestIntensityFraction(double[][] ionMatrix, int precursorCharge, Map<Double, Double> plMap, double ms2Tolerance) {
        int matchedPeakNum = 0;
        int maxRow = Math.max(2, Math.min(ionMatrix.length, 2 * (precursorCharge - 1)));
        int totalIonNum = ionMatrix[0].length * maxRow;
        Double[] intensityArray = plMap.values().toArray(new Double[0]);
        Arrays.sort(intensityArray, Collections.reverseOrder());
        double intensityT = 0;
        if (totalIonNum < intensityArray.length) {
            intensityT = intensityArray[totalIonNum];
        }
        int matchedHighestPeakNum = 0;
        for (int i = 0; i < maxRow; ++i) {
            for (int j = 0; j < ionMatrix[0].length; ++j) {
                for (double mz : plMap.keySet()) {
                    if (Math.abs(mz - ionMatrix[i][j]) <= ms2Tolerance) {
                        if (plMap.get(mz) > intensityT) {
                            ++matchedHighestPeakNum;
                        }
                        ++matchedPeakNum;
                        break;
                    }
                }
            }
        }

        if (matchedPeakNum > 0) {
            return (double) matchedHighestPeakNum / (double) matchedPeakNum;
        } else {
            return 0;
        }
    }

    public static double calExplainedAAFraction(double[][] ionMatrix, int precursorCharge, Map<Double, Double> plMap, double ms2Tolerance) {
        Set<Integer> matchedIdxSet = new HashSet<>();
        matchedIdxSet.add(0);
        matchedIdxSet.add(ionMatrix[0].length);
        int maxRow = Math.max(2, Math.min(ionMatrix.length, 2 * (precursorCharge - 1)));
        for (int i = 0; i < maxRow; ++i) {
            for (int j = 0; j < ionMatrix[0].length; ++j) {
                for (double mz : plMap.keySet()) {
                    if (Math.abs(mz - ionMatrix[i][j]) <= ms2Tolerance) {
                        if (i % 2 == 0) {
                            matchedIdxSet.add(j + 1);
                        } else {
                            if (j > 0) {
                                matchedIdxSet.add(j);
                            }
                        }
                        break;
                    }
                }
            }
        }


        Integer[] matchedIdxArray = matchedIdxSet.toArray(new Integer[0]);
        Arrays.sort(matchedIdxArray);
        int explainedAaNum = 0;
        if (matchedIdxArray.length > 1) {
            for (int i = 0; i < matchedIdxArray.length - 1; ++i) {
                if (matchedIdxArray[i + 1] - matchedIdxArray[i] == 1) {
                    ++explainedAaNum;
                }
            }
        }
        return (double) explainedAaNum / (double) ionMatrix[0].length;
    }

    public static double calAScore(TreeMap<Double, Double> plMap, int topN, Binomial binomial, TreeMap<Coordinate, Double> varPTMMap1, double[][] ionMatrix1, TreeMap<Coordinate, Double> varPTMMap2, double[][] ionMatrix2, double ms2Tolerance, int peptideLengthWithNC) throws Exception {
        double finalAScore = -9999;
        for (int localTopN = 1; localTopN <= topN; ++localTopN) {
            TreeMap<Double, Double> localPLMap = PrepareSpectrum.topNStyleNormalization(plMap, localTopN);
            double aScore = calAScoreSub(localPLMap, localTopN, binomial, varPTMMap1, ionMatrix1, varPTMMap2, ionMatrix2, ms2Tolerance, peptideLengthWithNC);
            if (aScore > finalAScore) {
                finalAScore = aScore;
            }
        }
       return finalAScore;
    }

    public static double calAScoreSub(TreeMap<Double, Double> localPLMap, int localTopN, Binomial binomial, TreeMap<Coordinate, Double> varPTMMap1, double[][] ionMatrix1, TreeMap<Coordinate, Double> varPTMMap2, double[][] ionMatrix2, double ms2Tolerance, int peptideLengthWithNC) throws Exception {
        TreeSet<Integer> totalAffectedBIonSet = new TreeSet<>();
        TreeSet<Integer> totalAffectedYIonSet = new TreeSet<>();
        if (varPTMMap2 == null) {
            getAffectedIonSet(varPTMMap1, peptideLengthWithNC - 2, totalAffectedBIonSet, totalAffectedYIonSet);
            Set<String> topMatchedPeakSet = getMatchedIonSet(ionMatrix1, localPLMap, ms2Tolerance, totalAffectedBIonSet, totalAffectedYIonSet);
            return  -10 * Math.log10(binomial.calProbLargerThanOrEqualTo(totalAffectedBIonSet.size() + totalAffectedYIonSet.size(), topMatchedPeakSet.size(), localTopN * 0.01));
        } else {
            getAffectedIonSet(varPTMMap1, peptideLengthWithNC - 2, totalAffectedBIonSet, totalAffectedYIonSet);
            getAffectedIonSet(varPTMMap2, peptideLengthWithNC - 2, totalAffectedBIonSet, totalAffectedYIonSet);


            if (!totalAffectedBIonSet.contains(1)) {
                totalAffectedBIonSet.pollFirst();
            }
            if (!totalAffectedBIonSet.contains(peptideLengthWithNC - 3)) {
                totalAffectedBIonSet.pollLast();
            }
            if (!totalAffectedYIonSet.contains(1)) {
                totalAffectedYIonSet.pollFirst();
            }
            if (!totalAffectedYIonSet.contains(peptideLengthWithNC - 3)) {
                totalAffectedYIonSet.pollLast();
            }

            Set<String> topMatchedPeakSet = getMatchedIonSet(ionMatrix1, localPLMap, ms2Tolerance, totalAffectedBIonSet, totalAffectedYIonSet);
            Set<String> secondMatchedPeakSet = getMatchedIonSet(ionMatrix2, localPLMap, ms2Tolerance, totalAffectedBIonSet, totalAffectedYIonSet);
            return -10 * Math.log10(binomial.calProbLargerThanOrEqualTo(totalAffectedBIonSet.size() + totalAffectedYIonSet.size(), topMatchedPeakSet.size(), localTopN * 0.01)) + 10 * Math.log10(binomial.calProbLargerThanOrEqualTo(totalAffectedBIonSet.size() + totalAffectedYIonSet.size(), secondMatchedPeakSet.size(), localTopN * 0.01));
        }
    }

    public static int getMatchedIonNum(TreeMap<Double, Double> plMap, int localMaxMs2Charge, double[][] ionMatrix, double ms2Tolerance) {
        int K1 = 0;
        for (int i = 0; i < localMaxMs2Charge * 2; ++i) {
            for (int j = 0; j < ionMatrix[0].length; ++j) {
                for (double mz : plMap.keySet()) {
                    if (Math.abs(mz - ionMatrix[i][j]) <= ms2Tolerance) {
                        ++K1;
                        break;
                    }
                }
            }
        }
        return K1;
    }

    private static void getAffectedIonSet(TreeMap<Coordinate, Double> varPtmMap, int peptideLength, Set<Integer> affectedBIonSet, Set<Integer> affectedYIonSet) {
        for (Coordinate co : varPtmMap.keySet()) {
            if (co.x == 0 || co.x == 1) {
                affectedBIonSet.add(1);
                affectedYIonSet.add(peptideLength - 1);
            } else if (co.x == peptideLength + 1 || co.x == peptideLength) {
                affectedBIonSet.add(peptideLength - 1);
                affectedYIonSet.add(1);
            } else {
                affectedBIonSet.add( co.x - 1);
                affectedBIonSet.add(co.x);
                affectedYIonSet.add( peptideLength - co.x);
                affectedYIonSet.add(peptideLength - co.x + 1);
            }
        }
    }

    private static Set<String> getMatchedIonSet(double[][] ionMatrix, TreeMap<Double, Double> localPLMap, double ms2Tolerance, Set<Integer> affectedBIonSet, Set<Integer> affectedYIonSet) {
        Set<String> matchedIonSet = new HashSet<>();
        for (int ion : affectedBIonSet) {
            double theoMass = ionMatrix[0][ion - 1];
            for (double mz : localPLMap.keySet()) {
                if (Math.abs(mz - theoMass) <= ms2Tolerance) {
                    matchedIonSet.add(String.format(Locale.US, "b%d", ion));
                    break;
                }
            }
        }

        for (int ion : affectedYIonSet) {
            double theoMass = ionMatrix[1][ionMatrix[0].length - ion];
            for (double mz : localPLMap.keySet()) {
                if (Math.abs(mz - theoMass) <= ms2Tolerance) {
                    matchedIonSet.add(String.format(Locale.US, "y%d", ion));
                    break;
                }
            }
        }

        return matchedIonSet;
    }
}
