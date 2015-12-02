import java.lang.Math;
import java.util.ArrayList;
//NEED TO RECONFIGURE//
public class Distance{
  //METHODS///////////////////////////////
  ////euclideanDistance(Atom a1, Atom a2);
  ////manhattanDistance(Coords c1);
  ////distanceSquared(Coords c1);
  ////lcs(); *longest common sequence
  ////lrmsd(Model1, Model2);
  ////lrmsd(Structure n[ative], Structure confs);
  ////////////////////////////////////////
  public static double euclideanDistance(Atom a1, Atom a2){
    return Math.sqrt(distanceSquared(a1, a2));
  }
  
  public static double manhattanDistance(double[] c1, double[] c2){
    double deltaX = c1[0] - c2[0];
    double deltaY = c1[1] - c2[1];
    double deltaZ = c1[2] - c2[2];
    double result = Math.abs(deltaX) + Math.abs(deltaY) + Math.abs(deltaZ);
    return result;
  }
  //Distance between two atoms
  public static double distanceSquared(Atom a1, Atom a2){
    double[] c1 = a1.getCoords();
    double[] c2 = a2.getCoords();
    double deltaX = c1[0] - c2[0];
    double deltaY = c1[1] - c2[1];
    double deltaZ = c1[2] - c2[2];
    double result = (deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ);
    return result;
  }
  //Distance between two models
   public static double distanceSquared(Model m1, Model m2){
     ArrayList<Atom> a1 = m1.getAtomsList();
     ArrayList<Atom> a2 = m2.getAtomsList();
     double distance = 0.0;
     if(a1.size() == a2.size()){
       for(int i=0; i<a2.size(); i++){
         distance += distanceSquared(a1.get(i), a2.get(i));
       }
     }
     return distance/a1.size();
  }
   
   //public static Model[] lcs(Model m1, Model m2){
     //find longest common sequence
      //must also store direction - which one of the cases
      //in the recursive formulation gave you the LCS
      //CODE FROM AMARDA
     //create charater arrays of atom names
//     Model[] reducedModels = new Model[2];
//     char[] x = m1.getAtomNames();
//     char[] y = m2.getAtomNames();
//     //create matrix C of m+1 rows and n+1 cols, where
//     //m is the number of chars in sequence x and n in sequence y
//     char[][] C = new char[x.length+1][y.length+1];
//      if(x[i] == y[j]){
//        C[i, j] = C[i-1, j-1] + 1;
//        direction[i,j] = 1; //diagonal
//      }else{
//        C[i,j] = max{ C[i-1, j], C[i, j-1] }
//        if(C[i,j] = C[i-1, j]){
//          direction[j,j] = 0; //up
//        }else{ 
//          direction[i.j] = -1; //down
//          i = m;
//          j = n;
//        }
//      while(i > 0 && j > 0){
//        if(direction[i,j] == 1){
//          LCS[C[i,j]] = x[i];
//          positionfromxforcomparison = i;
//          positionfromyforcomparison = j;
//           }
//         else if(direction[i,j] == 0){
//          LCS[C[i,j]] = x[i]
//         }
//         else{
//          LCS[C[i,j]] = y[j]   //check if this is x[i] and above y[j]
//         }
//      }
    //}
      //return new Models with updated atoms lists, according to lcs
      //return reducedModels;
   //}
  
  ////lrmsd(Model m1);
  public static double lrmsd(Model m1, Model m2){
    double lrmsd;
    //realign to centroids
    Atom[] conf1_prime = realign(m1);
    Atom[] conf2_prime = realign(m2);
    
    //find best rotation of m1 with respect to m2
    //TBD
    
    //sqrt(abs(distance with rotation)/n)
    lrmsd = Math.sqrt(Math.abs(distanceSquared(m1, m2)/m1.size()));
    return lrmsd;
  }
  
  ////lrmsd(Structure n[ative], confs);
  public static ArrayList<Double> lrmsd(Structure n, Structure confs){
    if(n.getModel(0).size() != confs.getModel(0).size()){
        //compute list of positions for lcs then use to compute lrmsd
      }
    Model nativeModel = n.getModel(0);
    double currLRMSD = 0.0;
    ArrayList<Double> allRMSDs = new ArrayList<Double>();
    for(int i=0; i<confs.size(); i++){
      currLRMSD = lrmsd(nativeModel, confs.getModel(i));
      allRMSDs.add(currLRMSD);
    }
    return allRMSDs;
  }
  
  ////centroids(); **Just Alpha Carbons**
  public static double[] centroid(Model m1){
    double[] center = {0,0,0};
    double length = m1.getAlphaCarbons().length;
    for(int i=0; i < length; i++){
      center[0] += m1.getAlphaCarbons()[i].getCoords()[0];
      center[1] += m1.getAlphaCarbons()[i].getCoords()[1];
      center[2] += m1.getAlphaCarbons()[i].getCoords()[2];
    }
    for(int j=0; j<3; j++){
      center[j] /= length;
    }
    return center;
  }
  ////realign();
  //Get each Atom's coords, subtract centroid, update coords
  public static Atom[] realign (Model m1){
    double avgX = centroid(m1)[0];
    double avgY = centroid(m1)[1];
    double avgZ = centroid(m1)[2];
    Atom[] caRealigned = m1.getAlphaCarbons();
    for(int i=0; i<caRealigned.length; i++){
      double x = caRealigned[i].getCoords()[0];
      double y = caRealigned[i].getCoords()[1];
      double z = caRealigned[i].getCoords()[2];
      caRealigned[i].setCoords((x-avgX), (y-avgY), (z-avgZ));
    }
    return caRealigned;
  } 
}
