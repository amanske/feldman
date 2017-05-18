import java.math.*;
import java.util.*;
import java.util.Scanner;
import java.io.*;

// Implemented using theory from https://en.wikipedia.org/wiki/Verifiable_secret_sharing, 
//https://profs.info.uaic.ro/~siftene/Feldman.pdf and https://courses.csail.mit.edu/6.857/2008/handouts/2005-Lecture16.pdf

public class Feldman {

	public static void main(String[] args) {
		BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
		String line;
		try {
			while ((line = br.readLine()) != null) {
				/////////////////// Read and store input //////////////////

				String[] data = line.split(" ");
				BigInteger p = new BigInteger(data[0]);
				BigInteger g = new BigInteger(data[1]);
				int d = Integer.parseInt(data[2]);
				int index = 4+d; // used to fetch from the right index in the data array
				ArrayList<BigInteger> a_i = new ArrayList<BigInteger>();

				for (int i = 3; i < index; i++) { // Start from index 3 (placement of d) and add all elements up to 3+d to a_i
					a_i.add(new BigInteger(data[i]));
				}

				int k = Integer.parseInt(data[index]); //get k (which is placed at index 4+d)
				index++; //increment by 1 to start iterating through k values

				ArrayList<BigInteger> s_i = new ArrayList<BigInteger>();
				for (int i = index; i < (index + k) ;i++) { // adds k elements to s_i
					s_i.add(new BigInteger(data[i]));
				}

				////////////////////// Parsing done //////////////////////

				// Now we have to find what s'_i = s_i
				// correctness is verified if g^s'_j mod p EQUALS bigPi_{i=0}^{d} (A_i^{j^i}) mod p

				ArrayList<BigInteger> ycoord = new ArrayList<BigInteger>();		// s_j
				ArrayList<Integer> xcoord = new ArrayList<Integer>();			// j

				
				for (int i = 1; i <= k; i++) { // iterate through indexes of s [1...k]
					BigInteger currentS = s_i.get(i-1); //get s_i
					BigInteger leftEq = g.modPow(currentS, p);
					BigInteger rightEq = calculateSum(a_i, i, p, d);
					if (leftEq.compareTo(rightEq) == 0) {	// correctness verified
						ycoord.add(currentS); //add s_i to y coordinate
						xcoord.add(i); // add i to x coordinate
					}
					if (xcoord.size() >= d+1) break; 	// enough points found, break loop
				}

				// Coordinates aquired, calculate a_0 with lagrange interpolation
				// Theory from MIT ppt
				BigInteger result = BigInteger.ZERO;
				BigInteger q = BigInteger.valueOf((p.intValue()-1)/2); 
				for (int i = 0; i <= d; i++) {
					BigInteger lower = BigInteger.ONE;
					BigInteger upper = BigInteger.ONE;
					BigInteger term = ycoord.get(i);
					for (int j = 0; j <= d; j++) {
						if (i != j) {
							upper = upper.multiply(BigInteger.valueOf(xcoord.get(j)).negate()).mod(q);
							BigInteger lowerpart = BigInteger.valueOf(xcoord.get(i)).subtract(BigInteger.valueOf(xcoord.get(j)));
							lower = lower.multiply(lowerpart).mod(q);
						}
					}
					BigInteger inverseLower = lower.modInverse(q);
					term = term.multiply(inverseLower).multiply(upper).mod(q);
					result = result.add(term).mod(q);
				}

				System.out.println(result.intValue()); //Print result

			}
		} catch (IOException e) {

		}
	}

	// Caluculate the right hand sum in the verifying equation
	public static BigInteger calculateSum(ArrayList<BigInteger> a_i, int j, BigInteger prime, int d) {
		BigInteger sum = BigInteger.ONE;
		BigInteger nextValue;
		BigInteger jindex = new BigInteger(Integer.toString(j));
		for (int i = 0; i <= d; i++) {
			BigInteger exponent = jindex.pow(i); 
			nextValue = a_i.get(i).modPow(exponent, prime);	//a_i ^ j^ i % p
			sum = sum.multiply(nextValue).mod(prime);
		}
		return sum;
	}
}