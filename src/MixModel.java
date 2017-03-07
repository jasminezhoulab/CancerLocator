import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.BlockRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.FastMath;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class MixModel {
	private RealMatrix[] mixDens;
	private int nBetas; // # of pre-calculated beta values
	private int nPoints=201; // # of points in integral approximation
	public MixModel(MethyModel tumor, MethyModel normal, RealVector thetas, int nBetas, int MYTHREADS) throws InterruptedException {
		int nFeatures=tumor.getNaRatio().getDimension();
		this.nBetas = nBetas;
		RealVector betas = new ArrayRealVector(nBetas);
		for (int i=0; i<nBetas; i++) {
			betas.setEntry(i,i/(nBetas-1.0));
		}
		mixDens = new RealMatrix[nFeatures];
		ExecutorService executor = Executors.newFixedThreadPool(MYTHREADS);

		for(int i = 0; i < nFeatures; i++) {
			double tumorAlpha = tumor.getAlpha().getEntry(i);
			double tumorBeta = tumor.getBeta().getEntry(i);
			BetaDistribution tumorDist = new BetaDistribution(tumorAlpha,tumorBeta);
			double normalAlpha = normal.getAlpha().getEntry(i);
			double normalBeta = normal.getBeta().getEntry(i);
			BetaDistribution normalDist = new BetaDistribution(normalAlpha,normalBeta);
			Runnable worker = new CalMixDens(tumorDist,normalDist,thetas,betas,nPoints,i,mixDens);
			executor.execute(worker);
		}
		executor.shutdown();
		while (!executor.isTerminated()) {
			Thread.sleep(10000);
		}
	}
    
	public MixModel() {		
	}
	
	public MixModel selectFeature(boolean[] selectedFeatures) {
		MixModel newModel = new MixModel();
		int nSelectedFeature = 0;
		for (boolean select:selectedFeatures) {
			if (select) {nSelectedFeature++;}
		}
		RealMatrix[] newMixDens = new RealMatrix[nSelectedFeature];
		int j = 0;
		for(int i=0; i<selectedFeatures.length; i++) {
			if (selectedFeatures[i]) {
				newMixDens[j] = this.mixDens[i];
				j++;
			}
		}
		newModel.setMixDens(newMixDens);
		return newModel;
	}
	
	public MixModel selectFeature(int[] selectedFeatures) {
		MixModel newModel = new MixModel();
		int nSelectedFeature = selectedFeatures.length;
		RealMatrix[] newMixDens = new RealMatrix[nSelectedFeature];
		for(int i=0; i<nSelectedFeature; i++) {
				newMixDens[i] = this.mixDens[selectedFeatures[i]];
		}
		newModel.setMixDens(newMixDens);
		return newModel;
	}

	public static class CalMixDens implements Runnable {
		BetaDistribution tumorDist;
		BetaDistribution normalDist;
		RealVector thetas;
		RealVector betas;
		int nThetas;
		int nBetas;
		int nPoints;
		int featureIdx;
		RealMatrix[] mixDens;

		CalMixDens(BetaDistribution tumorDist, BetaDistribution normalDist, RealVector thetas, RealVector betas,
				   int nPoints, int featureIdx, RealMatrix[] mixDens) {
			this.tumorDist = tumorDist;
			this.normalDist = normalDist;
			this.thetas = thetas;
			this.nThetas = thetas.getDimension();
			this.betas = betas;
			this.nBetas = betas.getDimension();
			this.nPoints = nPoints;
			this.featureIdx = featureIdx;
			this.mixDens = mixDens;
		}

		@Override
		public void run() {
			mixDens[featureIdx] = new BlockRealMatrix(nBetas,nThetas);
			for(int j = 0; j < nThetas;j++) {
				double theta = thetas.getEntry(j);
				RealVector betaDens = new ArrayRealVector(nBetas);
				for (int k = 0; k < nBetas; k++) {
					double beta = betas.getEntry(k);
					double lowerBound = FastMath.max(0, (beta - 1 + theta) / theta);
					if (Double.isNaN(lowerBound)) lowerBound = 0;
					double upperBound = FastMath.min(1,beta/theta);
					if (Double.isNaN(upperBound)) upperBound = 1;
					double step = (upperBound-lowerBound)/(nPoints-1);
					RealVector dens = new ArrayRealVector(nPoints);
					RealVector points = new ArrayRealVector(nPoints);
					RealVector allTumorDens = new ArrayRealVector(nPoints);
					RealVector allNormalDens = new ArrayRealVector(nPoints);
					RealVector allNormalDensRev = new ArrayRealVector(nPoints);
					// tumor
					for (int l = 0; l < nPoints; l++) {
						double tumorValue = lowerBound + l * step;
						points.setEntry(l, tumorValue);
						if (tumorValue == 0) tumorValue = 0.0001;
						if (tumorValue == 1) tumorValue = 0.9999;
						allTumorDens.setEntry(l,tumorDist.density(tumorValue));
					}
					// adjust the densities
					double calProb = tumorDist.probability(lowerBound,upperBound);
					double estProb = CancerLocator.integSimpson(points,allTumorDens);
					if (estProb!=0)	{
						allTumorDens.mapMultiplyToSelf(calProb/estProb);
					}
					else {
						allTumorDens.mapAddToSelf(1.0/allTumorDens.getDimension());
					}
					// normal
					RealVector normalPoints = new ArrayRealVector(nPoints);
					for (int l = 0; l < nPoints; l++) {
						double normalValue = (beta-theta*points.getEntry(l))/(1-theta);
						normalPoints.setEntry(nPoints-l-1,normalValue);
						if (normalValue == 0) normalValue = 0.0001;
						if (normalValue == 1) normalValue = 0.9999;
						double normalDens = normalDist.density(normalValue);
						allNormalDens.setEntry(l,normalDens);
						allNormalDensRev.setEntry(nPoints-l-1,normalDens);
					}
					calProb = normalDist.probability((beta-theta*upperBound)/(1-theta),(beta-theta*lowerBound)/(1-theta));
					estProb = CancerLocator.integSimpson(normalPoints, allNormalDensRev);
					if (estProb!=0)	{
						allNormalDens.mapMultiplyToSelf(calProb/estProb);
					}
					else {
						allNormalDens.mapAddToSelf(1.0/allNormalDens.getDimension());
					}
					//mixture
					for (int l = 0; l < nPoints; l++) {
						dens.setEntry(l,allTumorDens.getEntry(l)*allNormalDens.getEntry(l));
					}
					betaDens.setEntry(k,CancerLocator.integSimpson(points,dens));
				}
				double normTerm = CancerLocator.integSimpson(betas,betaDens); //normalization term
				if (normTerm!=0) {
					betaDens.mapDivideToSelf(normTerm);
				}
				else {
					betaDens.mapAddToSelf(1.0/betaDens.getDimension());
				}
				mixDens[featureIdx].setColumnVector(j, betaDens);
			}
		}
	}

	public RealMatrix[] getMixDens() {
		return mixDens;
	}	
	
	public void setMixDens(RealMatrix[] mixDens) {
		this.mixDens = mixDens;
	}

	public int getFeatureNum() {
		return mixDens.length;
	}
	
	public int getThetaNum() {
		return mixDens[0].getColumnDimension();
	}


	public int getnBetas() {
		return nBetas;
	}

	public void setnBetas(int nBetas) {
		this.nBetas = nBetas;
	}
}
