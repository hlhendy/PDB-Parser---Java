import weka.core.Instances;
import java.io.BufferedReader;
import java.io.FileReader;
import weka.filters.supervised.instance.Resample;
import weka.filters.supervised.instance.SMOTE;
import weka.filters.supervised.instance.SpreadSubsample;

public static void main(String[] args){
	/////////////
	//Read File//
	/////////////

	//read file name - MUST BE .arff format
	BufferedReader reader = new BufferedReader(new FileReader(filename));
	Instances data = new Instances(reader);
	reader.close();
	data.setClassIndex(data.numAttributes()-1);
	/////////////////
	//Preprocessing//
	/////////////////

	//Resample: random subsample
	//OPTIONS
	//biasToUniformClass -- Whether to use bias towards a uniform class. A value of 0 leaves the class distribution as-is, a value of 1 ensures the class distribution is uniform in the output data.
	//invertSelection -- Inverts the selection (only if instances are drawn WITHOUT replacement).
	//noReplacement -- Disables the replacement of instances.
	//randomSeed -- Sets the random number seed for subsampling.
	//sampleSizePercent -- The subsample size as a percentage of the original set.
	Instances resampled = weka.filters.supervised.instance.Resample(0.0, False, False, 1, 100.0)

	//SMOTE: synthetic minority oversampling technique
	//OPTIONS
	//classValue -- The index of the class value to which SMOTE should be applied. Use a value of 0 to auto-detect the non-empty minority class.
	//nearestNeighbors -- The number of nearest neighbors to use.
	//percentage -- The percentage of SMOTE instances to create.
	//randomSeed -- The seed used for random sampling.
	Instances smoted = weka.filters.supervised.instance.SMOTE(0, 5, 100.0, 1)

	//SpreadSubsample
	//OPTIONS
	//adjustWeights -- Whether instance weights will be adjusted to maintain total weight per class.
	//distributionSpread -- The maximum class distribution spread. (
		   //0 = no maximum spread, 1 = uniform distribution, 10 = allow at most a 10:1 ratio between the classes).
	//maxCount -- The maximum count for any class value (0 = unlimited).
	//randomSeed -- Sets the random number seed for subsampling.
	Instances spreadsubsampled = weka.filters.supervised.instance.SpreadSubsample(False, 1.0, 0.0, 1)

	///////////////////////
	//Attribute Selection//
	///////////////////////

	//InfoGrain with Ranker; 10-fold cross-validation
	weka.attributeSelection.InfoGainAttributeEval
}
