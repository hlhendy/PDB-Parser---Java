//NEED TO RECONFIGURE//
public class Distance{
  //METHODS///////////////////////////////
  ////euclideanDistance(Coords c1);
  ////manhattanDistance(Coords c1);
  ////distanceSquared(Coords c1);
  ////////////////////////////////////////
  public static double euclideanDistance(double[] c1, double[] c2){
    return Math.sqrt(distanceSquared(c1, c2));
  }
  
  public static double manhattanDistance(double[] c1, double[] c2){
    double deltaX = c1[0] - c2[0];
    double deltaY = c1[1] - c2[1];
    double deltaZ = c1[2] - c2[2];
    double result = Math.abs(deltaX) + Math.abs(deltaY) + Math.abs(deltaZ);
    return result;
  }
  public static double distanceSquared(double[] c1, double[] c2){
    double deltaX = c1[0] - c2[0];
    double deltaY = c1[1] - c2[1];
    double deltaZ = c1[2] - c2[2];
    double result = (deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ);
    return result;
  }
  ////lrmsd(Model m1);
  public double lrmsd(Model m1, Model m2){
    Atom[] conf1_prime = realign(m1);
    Atom[] conf2_prime = realign(m2);
    
    //calculate optimal rotation
    return 0.0;
  }
  
  ////centroids(); **Just Alpha Carbons**
  public double[] centroid(Model m1){
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
  public Atom[] realign (Model m1){
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
  
  ////euclideanDistance(m1); *Array of atoms*
  public static double euclideanDistance(Model m1, Model m2){
    double distance = 0.0;
    Atom[] m1CA = m1.getAlphaCarbons();
    Atom[] m2CA = m2.getAlphaCarbons();
    for(int i=0; i < m1.getAlphaCarbons().length; i++){
      distance = distance + Distance.euclideanDistance(m1CA[i].getCoords(), m2CA[i].getCoords());
    }
    return distance;
  }
  
}