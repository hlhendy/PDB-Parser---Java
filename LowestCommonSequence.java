public class LowestCommonSequence{
  public static void main(String[] args){
    if(args.length != 2){
      System.out.println("Please enter two character arrays.");
      return;
    }
    char[] x = args[0].toCharArray();
    char[] y = args[1].toCharArray();
    int m = x.length;
    int n = y.length;
    int[][] C = new int[m+1][n+1];
    int[][] direction = new int[m+1][n+1];
    C[0][0] = 0;
    
    for(int i=m-1; i >=0 ; i--){
      for(int j=n-1; j>=0; j--){
        if((x[i] != '\0') && (y[j] != '\0')){
          if(x[i] == y[j]){
            C[i][j] = C[i+1][j+1] + 1;
            direction[i][j] = 1;
            System.out.println("Match");
          }
          else{
            C[i][j] = Math.max(C[i+1][j], C[i][j+1]);
          }
        }
      }
    }
    //recover LCS - modify with Amarda's text
    int i=0, j=0;
    while(i < m && j < n){
      if((x[i] != '\0') && (y[j] != '\0')){
        if(x[i] == y[j]){
          //do something
          i++;
          j++;
        }
      }
      else if(C[i+1][j] >= C[i][j+1]){ i++; }
      else { j++; }
    }
  }
}