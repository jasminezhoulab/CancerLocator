import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.util.CombinatoricsUtils;
import org.apache.commons.math3.util.FastMath;

import java.io.*;
import java.text.DecimalFormat;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.Stream;

public class CancerLocator {
	public static void main(String[] args) throws IOException, InterruptedException {
		//String currentDir = System.getProperty("user.dir");

		// command line args
		String configFile = args[0]; //  config file

		// read the config file
		Properties prop = new Properties();
		InputStream config;
		config = new FileInputStream(configFile);
		prop.load(config);
		// training samples
		String trainFile = prop.getProperty("trainFile");
		// testing samples, methylation value
		String testMethyFile = prop.getProperty("testMethyFile");
		// testing samples, reads depth
		String testDepthFile = prop.getProperty("testDepthFile");
		// sample type to tissue type mapping file
		String typeMappingFile = prop.getProperty("typeMappingFile");
		// prediction results
		String resultFile = prop.getProperty("resultFile");
		PrintWriter results = new PrintWriter(resultFile, "UTF-8");
		results.println(join(new String[] {"Sample ID", "Log-likelihood ratio", "Predicted tumor burden",
				"Predicted sample class"}, "\t"));
		// methylation range cutoff used in feature filtering
		double rangeCut = Double.parseDouble(prop.getProperty("methylationRangeCutoff"));
		// likelihood ratio cutoff used in prediction
		double ratioCut = Double.parseDouble(prop.getProperty("logLikelihoodRatioCutoff"));
		// theta step
		//double thetaStep = Double.parseDouble(prop.getProperty("thetaStep"));
		double thetaStep = 0.01;
		// #threads used
		int nThreads = Integer.parseInt(prop.getProperty("nThreads"));

		int nBetas = 201; // num of beta values used in Simpson integration

		System.out.println("Run configration:");
		for (String key : prop.stringPropertyNames()) {
			String value = prop.getProperty(key);
			System.out.println(key + ": " + value);
		}
		System.out.println();

		// sample type to prediction class mapping
		HashMap<String, String> type2Class = new HashMap<>();
		String normalType="normal"; // the type name of normal samples
		BufferedReader br = new BufferedReader(new FileReader(typeMappingFile));
		String line;
		while((line=br.readLine())!=null) {
			String[] fields = line.trim().split("\t");
			type2Class.put(fields[0],fields[1]);
			//normal samples may have alternative label in training file
			if (fields[1].equals("normal")) normalType=fields[0];
		}

		// the theta values
		int nThetas= (int) (1/thetaStep + 1);
		List<Double> thetaList = DoubleStream.iterate(0, n -> n + thetaStep).limit(nThetas)
				.boxed().collect(Collectors.toList());
		if (thetaList.get(nThetas-1) >= 0.9999) { //theta must be less than 1
			thetaList.remove(nThetas-1);
			nThetas--;
		}
		RealVector thetas = new ArrayRealVector(thetaList.toArray(new Double[nThetas]));

		// read the training data
		int featureNum = 0;
		Map<String, ArrayList<double[]>> trainBeta = new HashMap<>();
		br = new BufferedReader(new FileReader(trainFile));
		while ((line = br.readLine()) != null) {
			String[] fields = line.trim().split("\t");
			if (featureNum == 0) featureNum = fields.length - 1;
			String type = fields[0]; // type of the sample
			double[] betaValue = str2double(Arrays.copyOfRange(fields, 1, fields.length));
			ArrayList<double[]> betaValues = trainBeta.get(type);
			if (betaValues == null) {
				betaValues = new ArrayList<>();
				trainBeta.put(type, betaValues);
			}
			betaValues.add(betaValue);
		}
		br.close();

		// build methylation models
		List<String> sampleTypes = new ArrayList<>(trainBeta.keySet());
		List<String> diseaseTypes = new ArrayList<>(sampleTypes);
		diseaseTypes.remove(normalType);
		Map<String, RealMatrix> rawValues = new HashMap<>();
		Map<String, MethyModel> models = new HashMap<>();
		for (String type : sampleTypes) {
			ArrayList<double[]> data = trainBeta.get(type); //beta-values
			int sampleSize = data.size();
			RealMatrix dataMatrix = new BlockRealMatrix(sampleSize, featureNum);
			for (int i = 0; i < sampleSize; i++) {
				dataMatrix.setRow(i, data.get(i));
			}
			models.put(type, new MethyModel(dataMatrix));
			rawValues.put(type, dataMatrix);
		}

		// filter features
		boolean[] selectedFeatures = new boolean[featureNum];
		Arrays.fill(selectedFeatures, Boolean.TRUE);

		// within each type
		for (String type : sampleTypes) {
			RealVector alpha = models.get(type).getAlpha();
			RealVector beta = models.get(type).getBeta();
			for (int i = 0; i < featureNum; i++) {
				if (!(alpha.getEntry(i) > 0 && beta.getEntry(i) > 0)) {
					selectedFeatures[i] = false;
				}
			}
		}

		// across types
		Map<Integer, Double> featureRanges = new HashMap<>();
		for (int i = 0; i < featureNum; i++) {
			if (selectedFeatures[i]) {
				double ncr = 0; // normal cancer range
				double ctr = 0; // cancer type range
				double[] means = new double[sampleTypes.size()];

				for (int typeInd = 0; typeInd < sampleTypes.size(); typeInd++) {
					String type = sampleTypes.get(typeInd);
					double[] values = rawValues.get(type).getColumn(i);
					means[typeInd] = calMean(values);
				}

				for (int typeInd1 = 0; typeInd1 < sampleTypes.size(); typeInd1++) {
					String type1 = sampleTypes.get(typeInd1);
					for (int typeInd2 = typeInd1 + 1; typeInd2 < sampleTypes.size(); typeInd2++) {
						String type2 = sampleTypes.get(typeInd2);
						double methyDiff = Math.abs(means[typeInd1]-means[typeInd2]);
						if (type1.equals(normalType) | type2.equals(normalType)) {
							if (methyDiff > ncr) {
								ncr = methyDiff;
							}
						} else {
							if (methyDiff > ctr) {
								ctr = methyDiff;
							}
						}
					}
				}

				double featureRange = Math.max(ncr, ctr);
				featureRanges.put(i, featureRange);
			}
		}

		int nSelectedFeatures = 0;
		for (int i = 0; i < featureNum; i++) {
			if (selectedFeatures[i]) {
				if (featureRanges.get(i) < rangeCut) {
					selectedFeatures[i] = false;
				}
				else {
					nSelectedFeatures++;
				}
			}
		}

		// build the mixture models
		// only selected features considered
		System.out.println("Calculating the mixture models...");
		Map<String, HashMap<Integer, MixModel>> mixModels = new HashMap<>();
		boolean[] goodMixModels = new boolean[nSelectedFeatures];
		Arrays.fill(goodMixModels, Boolean.TRUE);
		MethyModel ctrModel = models.get(normalType);
		for (String type : diseaseTypes) {
			MethyModel diseaseModel = models.get(type);
			HashMap<Integer, MixModel> typeMixModels = new HashMap<>();
			for (int copyNum : new int[]{2}) { // no CNV considered
				// thetas at DNA level
				RealVector thetasDNA = calThetasDNA(thetas, copyNum);
				MixModel mix = new MixModel(diseaseModel.selectFeature(selectedFeatures),
						ctrModel.selectFeature(selectedFeatures), thetasDNA, nBetas, nThreads);
				typeMixModels.put(copyNum, mix);
				for (int featureIdx = 0; featureIdx < nSelectedFeatures; featureIdx++) {
					RealMatrix mixDens = mix.getMixDens()[featureIdx];
					for (int beta = 0; beta < mixDens.getRowDimension(); beta++) {
						double dens = mixDens.getEntry(beta, 0);
						if (Double.isNaN(dens) || Double.isInfinite(dens)) {
							goodMixModels[featureIdx] = false;
						}
					}
				}
			}
			System.out.println(type+" model calculated");
			mixModels.put(type, typeMixModels);
		}
		System.out.println();

		// update mixture models
		for (String type : diseaseTypes) {
			MixModel model = mixModels.get(type).get(2).selectFeature(goodMixModels);
			mixModels.get(type).put(2, model);
		}

		// update good features
		int j = 0;
		for (int i = 0; i < featureNum; i++) {
			if (selectedFeatures[i]) {
				if (!goodMixModels[j]) {// don't use this feature
					selectedFeatures[i] = false;
					nSelectedFeatures--;
				}
				j++;
			}
		}
		System.out.println(nSelectedFeatures + " features used in inference");
		System.out.println();

		// load the testing data
		// read the depth file first
		Map<String, MethySample> testSet = new HashMap<>();
		br = new BufferedReader(new FileReader(testDepthFile));
		while ((line = br.readLine()) != null) {
			String[] fields = line.trim().split("\t");
			String sampleID = fields[0];
			String[] temp = Arrays.copyOfRange(fields, 1, fields.length);
			int[] depth = str2int(temp);
			MethySample sample = new MethySample(sampleID, fields[0], 0, -1, true);
			sample.setDepth(depth);
			testSet.put(sampleID, sample);
		}
		br.close();


		// then read the file with methylated CpG numbers
		br = new BufferedReader(new FileReader(testMethyFile));
		List<String> testSamples = new ArrayList<>();
		while ((line = br.readLine()) != null) {
			String[] fields = line.trim().split("\t");
			String sampleID = fields[0];
			testSamples.add(sampleID);
			String[] temp = Arrays.copyOfRange(fields, 1, fields.length);
			int[] methyDepth = str2int(temp);
			int[] depth = testSet.get(sampleID).getDepth();
			double[] methy = new double[methyDepth.length];
			for (int i = 0; i < methy.length; i++) {
				if (depth[i] > 0) {
					methy[i] = ((double) methyDepth[i]) / depth[i];
				} else {
					methy[i] = Double.NaN;
				}
			}
			testSet.get(sampleID).setMethy(methy);
			testSet.get(sampleID).selfSelectFeature(selectedFeatures);
		}
		br.close();

		System.out.println("Making predictions...");
		// make predictions on the test set
		Map<String, MethySample> predictions = new HashMap<>();
		//multithreading
		ExecutorService executor = Executors.newFixedThreadPool(nThreads);
		for (String sampleId : testSet.keySet()) {
			MethySample testSample = testSet.get(sampleId);
			Runnable worker = new Predictor(testSample, mixModels, predictions, diseaseTypes, thetas);
			executor.execute(worker);
		}
		executor.shutdown();

		while (!executor.isTerminated()) {
			Thread.sleep(10000);
		}

		// write the results
		// in the same order of samples
		for (String sampleId : testSamples) {
			MethySample predSample = predictions.get(sampleId);
			double predTheta;
			try {
				predTheta = predSample.getTheta();
			}  catch (Exception e) {
				predTheta = -1;
			}
			double densRatio = predSample.getDensRatio()/nSelectedFeatures; //normalized by feature number
			String predType = predSample.getType();
			if (densRatio==0) predType = normalType; // no matter what cutoff used
			String predClass = densRatio<ratioCut?type2Class.get(normalType):type2Class.get(predType);
			DecimalFormat thetaFormat = new DecimalFormat("#.###"); // for predicted theta
			String thetaOutput = thetaFormat.format(predTheta);
			String output = join(new String[] {sampleId, Double.toString(densRatio), thetaOutput, predClass} ,"\t");
			results.println(output);
		}
		results.close();
		System.out.println("All jobs done!");
	}

	private static double calMean(double[] values) {
		List<Double> nonNa = new ArrayList<>();
		for (int i = 0; i<values.length; i++) {
			if (!Double.isNaN(values[i])) {
				nonNa.add(values[i]);
			}
		}
		int nonNaCount = nonNa.size();
		double[] nonNaValues = Stream.of(nonNa.toArray(new Double[nonNaCount])).
				mapToDouble(Double::doubleValue).toArray();
		Mean mean = new Mean();
		return mean.evaluate(nonNaValues,0,nonNaValues.length);
	}

	private static MethySample samplePred(MethySample subTestSample,
										  Map<String, HashMap<Integer, MixModel>> subMixModels,
										  List<String> diseaseTypes, RealVector thetas) {

		String sampleId = subTestSample.getId();
		int nTheta = subMixModels.get(diseaseTypes.get(0)).get(2).getThetaNum();

		// prediction w/o CNA considered
		RealMatrix dens = new BlockRealMatrix(diseaseTypes.size(), nTheta);
		for (int i = 0; i < diseaseTypes.size(); i++) {
			String type = diseaseTypes.get(i);
			HashMap<Integer, MixModel> typeSubMixModels = subMixModels.get(type);
			RealMatrix thetaDists = calSampleDens(subTestSample, typeSubMixModels.get(2));
			dens.setRowVector(i, calSumLogDens(thetaDists));
		}
		CancerPrediction pred = new CancerPrediction(dens, thetas, false);

		Integer[] typeRanks = pred.getTypeRanks();
		String typePred = diseaseTypes.get(typeRanks[0]);
		MethySample predSample = new MethySample(typePred, pred.getBestTheta());
		predSample.setId(sampleId);
		predSample.setDensRatio(pred.getBestRatio());
        return predSample;
	}

	public static class Predictor implements Runnable {
		MethySample subTestSample;
		Map<String, HashMap<Integer, MixModel>> subMixModels;
		Map<String, MethySample> predictions;
		List<String> diseaseTypes;
		RealVector thetas;
		Predictor(MethySample subTestSample, Map<String, HashMap<Integer, MixModel>> subMixModels,
				  Map<String, MethySample> predictions, List<String> diseaseTypes, RealVector thetas) {
			this.subTestSample = subTestSample;
			this.subMixModels = subMixModels;
			this.predictions = predictions;
			this.diseaseTypes = diseaseTypes;
			this.thetas = thetas;
		}

		@Override
		public void run() {
			String sampleId = subTestSample.getId();
			MethySample predSample = samplePred(subTestSample, subMixModels, diseaseTypes, thetas);
			predictions.put(sampleId, predSample);
		}
	}

	// calculate DNA level thetas
	public static RealVector calThetasDNA(RealVector thetas, int copyNum) {
		RealVector thetasDNA = thetas.mapMultiply((double) copyNum);
		if (copyNum != 0) {
			RealVector ctrRatios = thetas.mapMultiply(-1);
			ctrRatios.mapAddToSelf(1);
			ctrRatios.mapMultiplyToSelf(2);
			RealVector temp = thetasDNA.add(ctrRatios);
			thetasDNA = thetasDNA.ebeDivide(temp);
		}
		return thetasDNA;
	}

	private static RealMatrix calSampleDens(MethySample sample, MixModel model) {
		double[] betaValues = sample.getMethy();
		int[] depths = sample.getDepth();
		RealMatrix[] mixDens = model.getMixDens();
		int nFeature = sample.getFeatureNum();
		int nTheta = model.getThetaNum();
		int nBetas = mixDens[0].getRowDimension();
		RealVector betas = new ArrayRealVector(nBetas);
		for (int i = 0; i < nBetas; i++)
			betas.setEntry(i,i/(nBetas-1.0));
		RealMatrix dens = new BlockRealMatrix(nFeature, nTheta);
		for (int i = 0; i < nFeature; i++)
			for (int j = 0; j < nTheta; j++) {
				RealVector betaDens = mixDens[i].getColumnVector(j);
				int methyCounts = (int) Math.round(depths[i]*betaValues[i]);
				dens.setEntry(i, j, calLogCompoundDens(betas, betaDens, depths[i],methyCounts));
			}
		return dens;
	}

	private static double calLogCompoundDens(RealVector betas, RealVector betaDens, int n, int k) {
		double logComb = CombinatoricsUtils.binomialCoefficientLog(n,k);
		int nBetas = betas.getDimension();
		RealVector dens = new ArrayRealVector(nBetas);
		for (int i = 0; i < nBetas; i++) {
			dens.setEntry(i, betaDens.getEntry(i) * FastMath.pow(betas.getEntry(i), k)
					* FastMath.pow(1 - betas.getEntry(i), n - k));
		}
		double prob = integSimpson(betas,dens);
		double logProb=(prob==0)?-1000:FastMath.log(prob); // avoid -Inf
		return logComb+logProb;
	}


	private static RealVector calSumLogDens(RealMatrix thetaDists) {
		int nFeature = thetaDists.getRowDimension();
		int nTheta = thetaDists.getColumnDimension();
		RealVector sumLogDens = new ArrayRealVector(nTheta);
		for (int i = 0; i < nTheta; i++) {
			double sum = 0;
			int nonNaFeature = 0;
			for (int j = 0; j < nFeature; j++) {
				Double logDen = thetaDists.getEntry(j, i);
				if (!Double.isNaN(logDen)) {
					nonNaFeature++;
					sum += logDen;
				}
			}
			if (nonNaFeature != 0) {
				sumLogDens.setEntry(i, sum/nonNaFeature);
			}
			else { //no feature available
				sumLogDens.setEntry(i, Double.NaN);
			}
		}
		return sumLogDens;
	}

	public static double integSimpson(RealVector points, RealVector dens) {
		double s;
		int n = points.getDimension() - 1; // # of intervals
		double h = points.getEntry(1) - points.getEntry(0); // the length of an
		// interval
		if (n == 2) {
			s = dens.getEntry(0) + 4 * dens.getEntry(1) + dens.getEntry(2);
		} else {
			s = dens.getEntry(0) + dens.getEntry(n);
			for (int i = 1; i < n; i += 2) {
				s += 2 * dens.getEntry(i);
			}
			for (int i = 2; i < n - 1; i += 2) {
				s += 4 * dens.getEntry(i);
			}
		}
		s = s * h / 3;
		return s;
	}

	public static int[] str2int(String[] strings) {
		int[] intArray = new int[strings.length];
		int i = 0;
		for (String str : strings) {
			str = str.trim();
			if (str.equals("NA")) {
				intArray[i] = 0; // just for the depth info
			} else {
				intArray[i] = Integer.parseInt(str.trim());
			}
			i++;
		}
		return intArray;
	}

	public static double[] str2double(String[] strings) {
		double[] doubleArray = new double[strings.length];
		int i = 0;
		for (String str : strings) {
			str = str.trim();
			if (str.equals("NA")) {
				doubleArray[i] = Double.NaN;
			} else {
				doubleArray[i] = Double.parseDouble(str);
			}
			i++;
		}
		return doubleArray;
	}

	public static String join(String[] strings, String delim) {
		StringBuilder sb = new StringBuilder();
		String loopDelim = "";
		for (String s : strings) {
			sb.append(loopDelim);
			sb.append(s);
			loopDelim = delim;
		}
		return sb.toString();
	}
}