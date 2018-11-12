package snapper;

import java.io.IOException;
import java.util.Arrays;

import org.junit.Test;

import beast.app.beauti.BeautiDoc;
import junit.framework.TestCase;
import snap.likelihood.COMPLEX;

public class TestMatrixExponentiatorCF extends TestCase {
	final static String TEST1_DIR = "test/snapper/test1/";
	final static String TEST2_DIR = "test/snapper/test2/";
	final static String TEST3_DIR = "test/snapper/test3/";
	

	COMPLEX [] readComplexArray(String fileName) throws IOException {
		String s = BeautiDoc.load(fileName);
		String [] strs = s.split("\n");
		COMPLEX [] array = new COMPLEX[strs.length];
		for (int i = 0; i < array.length; i++) {
			array[i] = parseComplex(strs[i]);
		}
		return array;
	}


	COMPLEX []readComplexMatrix(String fileName) throws IOException {
		String s = BeautiDoc.load(fileName);
		String [] strs = s.split("\n");
		COMPLEX [][] array = new COMPLEX[strs.length][];
		for (int i = 0; i < array.length; i++) {
			s = strs[i];
			String [] strs2 = s.split("\\)");
			array[i] = new COMPLEX[strs2.length];
			for (int j = 0; j < array[i].length; j++) {
				array[i][j] = parseComplex(strs2[j]);
			}
		}
		
		COMPLEX [] AA = new COMPLEX[array.length * array[0].length];
		for (int i = 0; i < array.length; i++) {			
			for (int j = 0; j < array[0].length; j++) {
				AA[i*array[0].length + j] = array[i][j];
			}
		}
		return AA;
	}
	
	double []readDoubleMatrix(String fileName) throws IOException {
		String s = BeautiDoc.load(fileName);
		String [] strs = s.split("\n");
		int n = strs.length;
		double [] array = new double[n*n];
		for (int i = 0; i < n; i++) {
			s = strs[i];
			String [] strs2 = s.trim().split("\\s+");
			for (int j = 0; j < n; j++) {
				array[i * n + j] = Double.parseDouble(strs2[j]);
			}
		}
		return array;
	}
	
	private COMPLEX parseComplex(String string) {
		int startRe = string.indexOf('(') + 1;
		int endRe = string.indexOf('e') + 3;
		while (string.charAt(endRe) != '+') {
			endRe++;
		}
		int startIm = endRe + 1;
		int endIm = string.indexOf('j');
		String re = string.substring(startRe,endRe);
		String im = string.substring(startIm,endIm);
		try {
			return new COMPLEX(Double.parseDouble(re), Double.parseDouble(im));
		} catch (NumberFormatException e) {
			return null;
		}
	}

	@Test
	public void testGetValuesTridag() throws IOException {
		COMPLEX [] AA = readComplexMatrix(TEST1_DIR + "AA.txt");

		COMPLEX [] A2c = readComplexArray(TEST1_DIR + "A2c.txt");
		COMPLEX [] a = readComplexArray(TEST1_DIR + "a.txt");
		COMPLEX [] b = readComplexArray(TEST1_DIR + "b.txt");
		COMPLEX [] c = readComplexArray(TEST1_DIR + "c.txt");
		COMPLEX [] d = readComplexArray(TEST1_DIR + "d.txt");
		
		MatrixExponentiator e = new MatrixExponentiator();
		COMPLEX [][] out = new COMPLEX[4][];
		out[0] = new COMPLEX[a.length]; init(out[0]);
		out[1] = new COMPLEX[b.length]; init(out[1]);
		out[2] = new COMPLEX[c.length]; init(out[2]);
		out[3] = new COMPLEX[d.length]; init(out[3]);
		
		e.getValuesTridag(AA, A2c, out);
		
		compare(out[0], a);
		compare(out[1], b);
		compare(out[2], c);
		compare(out[3], d);
	}

	@Test
	public void testGetValuesUBand() throws IOException {
		COMPLEX [] AA = readComplexMatrix(TEST1_DIR + "AA.txt");
		double [] LL = readDoubleMatrix(TEST1_DIR + "LL.txt");

		COMPLEX [] c_i = readComplexArray(TEST1_DIR + "c_i.txt");

		COMPLEX [] A2c = readComplexArray(TEST1_DIR + "A2c.txt");
		COMPLEX [] a = readComplexArray(TEST1_DIR + "a0.txt");
		COMPLEX [] b = readComplexArray(TEST1_DIR + "b0.txt");
		COMPLEX [] c = readComplexArray(TEST1_DIR + "c0.txt");
		COMPLEX [] d = readComplexArray(TEST1_DIR + "d0.txt");
		
		MatrixExponentiator e = new MatrixExponentiator();
		COMPLEX [][] out = new COMPLEX[4][];
		out[0] = new COMPLEX[a.length]; init(out[0]);
		out[1] = new COMPLEX[b.length]; init(out[1]);
		out[2] = new COMPLEX[c.length]; init(out[2]);
		out[3] = new COMPLEX[d.length]; init(out[3]);
		
		COMPLEX z = new COMPLEX(-699.8688082445777,1399.5917029301354);
	
		
		e.getValuesUBand(LL, AA, A2c, c_i, z, out);
		
		compare(out[0], a);
		compare(out[1], b);
		compare(out[2], c);
		compare(out[3], d);
	}

	@Test
	public void testUBandSolver() throws IOException {
		COMPLEX [] a = readComplexArray(TEST1_DIR + "a0.txt");
		COMPLEX [] b = readComplexArray(TEST1_DIR + "b0.txt");
		COMPLEX [] c = readComplexArray(TEST1_DIR + "c0.txt");
		COMPLEX [] d = readComplexArray(TEST1_DIR + "d0.txt");
		
		COMPLEX [] solvEven = readComplexArray(TEST1_DIR + "sol_even.txt");
		COMPLEX [] x = new COMPLEX[solvEven.length]; init(x);
		
		MatrixExponentiator e = new MatrixExponentiator();
		COMPLEX [][] in = new COMPLEX[4][];
		in[0] = a;
		in[1] = b;
		in[2] = c;
		in[3] = d;

		e.ubandSolver(in, x, 33);
		
		compare(x, solvEven);
	}
	
	@Test
	public void testTriDagSolver() throws IOException {
		COMPLEX [] a = readComplexArray(TEST1_DIR + "a.txt");
		COMPLEX [] b = readComplexArray(TEST1_DIR + "b.txt");
		COMPLEX [] c = readComplexArray(TEST1_DIR + "c.txt");
		COMPLEX [] d = readComplexArray(TEST1_DIR + "d.txt");
		
		COMPLEX [] solvOdd = readComplexArray(TEST1_DIR + "sol_odd.txt");
		COMPLEX [] x = new COMPLEX[solvOdd.length]; init(x);
		
		MatrixExponentiator e = new MatrixExponentiator();
		COMPLEX [][] in = new COMPLEX[4][];
		in[0] = a;
		in[1] = b;
		in[2] = c;
		in[3] = d;

		e.triDiagSolver(in, x, 33);
		
		compare(x, solvOdd);
	}

	@Test
	public void testFasterSolver() throws IOException {
		COMPLEX [] AA = readComplexMatrix(TEST1_DIR + "AA.txt");
		double [] LL = readDoubleMatrix(TEST1_DIR + "LL.txt");

		COMPLEX [] c_i = readComplexArray(TEST1_DIR + "c_i.txt");

		COMPLEX [] A2c = readComplexArray(TEST1_DIR + "A2c.txt");
		
		COMPLEX z = new COMPLEX(-699.8688082445777,1399.5917029301354);

		MatrixExponentiator e = new MatrixExponentiator();

		COMPLEX [] x = e.fasterSolver(LL, AA, A2c, c_i, z);
		
		COMPLEX [] solv = readComplexArray(TEST1_DIR + "solve.txt");
		compare(x, solv);
	}

	@Test
	public void testExpCF() throws IOException {
		int N = 33;

		MatrixExponentiator e = new MatrixExponentiator();
		e.initMem(N);
		double [] v = new double[N];
		Arrays.fill(v, 1.0);
		QMatrix Q = new QMatrix(N);
		double dt = 0.1;
		double u_ = 0.001;
		double v_ = 0.001;
		double [] a = new double[]{u_,-(u_ + v_),0,0};
		double [] b = new double[]{0,2,-2,0};
		Q.setQ(a, b);
		double [] Qexpected = readDoubleMatrix(TEST1_DIR + "Q.txt");
		compare(Qexpected, Q.Q);

		double [] solv = new double[] {6.27377023e+00, 1.01641829e+01, 7.83332873e+00, 4.67503450e+00, 2.44284918e+00, 9.87181265e-01, 3.49305166e-01, 9.53919313e-02, 2.28199240e-02, 4.20205496e-03, 6.78490321e-04, 8.40547742e-05, 9.14713971e-06, 7.60753774e-07, 5.57178962e-08, 3.10680883e-09, 1.58673691e-10, -8.78717848e-13, -8.64448521e-12, -6.42884437e-13, 4.84353951e-12, 3.03469205e-12, -1.52385764e-12, -5.03889921e-12, -6.27086328e-12, -4.78906219e-12, -2.96315858e-12, -9.96341618e-13, 2.05139472e-13, 3.59732682e-13, 2.76602908e-13, -7.37648180e-13, -1.30750063e-12};
		e.expCF(dt, Q.Q, v);
		compare(v, solv, 1e-7);
		
		Q.initA2Q(N);
		double [] A2 = readDoubleMatrix(TEST1_DIR + "A2.txt");
		compare(A2, Q.A2);
		
		double [] A2Q = readDoubleMatrix(TEST1_DIR + "A2Q.txt");
		compare(A2Q, Q.A2Q);

		
		COMPLEX [][] w = new COMPLEX[e.ci_real.length][v.length];
		for (int i = 0; i < v.length; i++) {
			for (int j = 0; j < e.ci_real.length; j++) {
				w[j][i] = new COMPLEX();
			}
		}
		
		
		
//	    w = np.zeros([len(c_real),len(v)],dtype=np.complex)
		Arrays.fill(v, 1.0);
		double [] LL = Q.Q;
				
		COMPLEX [] vi = new COMPLEX[v.length];
		for (int i = 0; i < vi.length; i++) {
			vi[i] = new COMPLEX(1/dt,0);
		}		
		
		COMPLEX [] AA = new COMPLEX[LL.length];
		// TODO: initialise AA
		for (int i = 0; i < AA.length; i++) {
			AA[i] = new COMPLEX();
		}		
		
		COMPLEX [] A2c = new COMPLEX[N];
//		QMatrix.initA2(N);
		for (int i = 0; i < vi.length; i++) {
			double sum = 0;
			for (int k = 0; k < N; k++) {
				sum += Q.A2[i*N + k];
			}
			A2c[i] = new COMPLEX(sum/dt, 0);
		}		
		
		COMPLEX [] z_iExp = readComplexArray(TEST2_DIR + "z_i.txt");
		
		for (int i = 0; i < MatrixExponentiator.ci_real.length; i++) {
			System.err.println("Round " + i);
			// vi = np.sum(w,axis=0)
			double zi_real = MatrixExponentiator.zi_real[i]/dt; 
			double zi_imag = MatrixExponentiator.zi_imag[i]/dt;
			
			assertEquals(z_iExp[i].m_fRe, zi_real, 2e-13);
			assertEquals(z_iExp[i].m_fIm, zi_imag, 2e-13);
			
			int x = 0;
			for (int j = 0; j < AA.length; j++) {
				//for (int k = 0; k < N; k++) {
					AA[x].m_fRe = Q.A2Q[x] - zi_real * Q.A2[x];
					AA[x].m_fIm =          - zi_imag * Q.A2[x];
					x++;
				//}
			}
			
			double [] LLExp = readDoubleMatrix(TEST2_DIR + "LL" + i +".txt");
			compare(LLExp, LL);
			COMPLEX [] viExp = readComplexArray(TEST2_DIR + "c_i" + i +".txt");
			compare(viExp, vi);
			COMPLEX [] A2cExp = readComplexArray(TEST2_DIR + "A2c" + i +".txt");
			compare(A2cExp, A2c);
			COMPLEX [] AAExp = readComplexMatrix(TEST2_DIR + "AA" + i +".txt");
			compare(AAExp, AA);
			
			COMPLEX z = new COMPLEX(zi_real, zi_imag);
			
			// e.fasterSolver(LLExp, AAExp, A2cExp, vi, z);
			e.getValuesUBand(LL, AA, A2c, vi, z, e.out1);
			
			COMPLEX [] array = readComplexArray(TEST2_DIR + "uband0" + i +".txt");
			compareShortest(array, e.out1[0]);
			array = readComplexArray(TEST2_DIR + "uband1" + i +".txt");
			compareShortest(array, e.out1[1]);
			array = readComplexArray(TEST2_DIR + "uband2" + i +".txt");
			compareShortest(array, e.out1[2]);
			array = readComplexArray(TEST2_DIR + "uband3" + i +".txt");
			compareShortest(array, e.out1[3]);
			
			e.ubandSolver(e.out1, e.solvEven, N);
			array = readComplexArray(TEST2_DIR + "sol_even" + i +".txt");
			compareShortest(array, e.solvEven);

			e.getValuesTridag(AA, A2c, e.out2);
			array = readComplexArray(TEST2_DIR + "triband0" + i +".txt");
			compareShortest(array, e.out2[0]);
			array = readComplexArray(TEST2_DIR + "triband1" + i +".txt");
			compareShortest(array, e.out2[1]);
			array = readComplexArray(TEST2_DIR + "triband2" + i +".txt");
			compareShortest(array, e.out2[2]);
			array = readComplexArray(TEST2_DIR + "triband3" + i +".txt");
			compareShortest(array, e.out2[3]);
			
			e.triDiagSolver(e.out2, e.solvOdd, N);
			array = readComplexArray(TEST2_DIR + "sol_odd" + i +".txt");
			compareShortest(array, e.solvOdd, 1e-6);

			for (int k = A2c.length - 1; k >= 1; k -= 2) {
				e.solvEven[k].m_fRe = e.solvEven[k/2].m_fRe;
				e.solvEven[k].m_fIm = e.solvEven[k/2].m_fIm;
				e.solvEven[k-1].m_fRe = e.solvOdd[k/2-1].m_fRe;
				e.solvEven[k-1].m_fIm = e.solvOdd[k/2-1].m_fIm;
			}

			COMPLEX [] xExp = readComplexArray(TEST2_DIR + "x" + i +".txt");
			compare(xExp, e.solvEven, 1e-8);

			
			for (int j = 0; j < N; j++) {
				w[i][j].mul(e.solvEven[j], e.ci_real[i], e.ci_imag[i]);
			}
			COMPLEX [] wExpected = readComplexMatrix(TEST2_DIR + "w" + i +".txt");
			compare(wExpected, w, 1e-7);
		}
		
		w = new COMPLEX[v.length][e.ci_real.length];
		for (int i = 0; i < v.length; i++) {
			for (int j = 0; j < e.ci_real.length; j++) {
				w[i][j] = new COMPLEX();
			}
		}
		e.expCF(v, Q, dt, w);
		//COMPLEX [] wExpected = readComplexMatrix(TEST1_DIR + "w.txt");
		//compare(wExpected, w, 1e-7);

		
		

		e.expCF(v, Q, dt);
		compare(v, solv, 1e-7);

	}
	
	@Test
	public void testExpCF3() throws IOException {
		int N = 33;

		MatrixExponentiator e = new MatrixExponentiator();
		e.initMem(N);
		double [] v = new double[N];
		Arrays.fill(v, 1.0);
		QMatrix Q = new QMatrix(N);
		double dt = 0.01;
		double u_ = 0.5;
		double v_ = 0.5;
		double [] a = new double[]{u_,-(u_ + v_),0,0};
		double [] b = new double[]{0,2,-2,0};
		Q.setQ(a, b);

		double [] solv = new double[] {6.46372099e-03, -2.53441817e-03, -1.12124939e-02, 6.39070763e-03, 7.07907284e-03, -7.49556373e-03, -2.71519609e-03, 6.12625130e-03, -1.66842601e-04, -3.74391581e-03, 1.22319602e-03, 1.70014038e-03, -1.11547485e-03, -5.21067457e-04, 6.51493218e-04, 5.56120317e-05, -2.79146776e-04, 4.73038355e-05, 8.93791482e-05, -3.68112586e-05, -2.02636837e-05, 1.55215264e-05, 2.43339381e-06, -4.64622800e-06, 3.22485500e-07, 1.04019478e-06, -2.69229540e-07, -1.73318450e-07, 8.58345327e-08, 2.05479822e-08, -2.07225604e-08, -1.35330441e-09, 2.58243417e-09};
		v = new double[]{0.006463720994063693, -0.0025598894967633185, -0.011670084506784329, 0.006992547947822397, 0.008307368929750345, -0.009624494344782747, -0.003891770417487691, 0.009999979361246672, -3.164138031256699E-4, -0.008415978246663592, 0.003324991505858777, 0.0057013946826780845, -0.004708080014436136, -0.0028239150287447754, 0.004625163434420616, 5.276322657217052E-4, -0.0036109917030115145, 8.51152557681401E-4, 0.002282182301589475, -0.0013607669299408828, -0.0011063596425523177, 0.0012769476496068285, 3.077497378592821E-4, -9.215487768382256E-4, 1.0234023463673997E-4, 5.388342382900541E-4, -2.3224876169316427E-4, -2.5401044059004865E-4, 2.1803729621996514E-4, 9.22966179980842E-5, -1.679166495819626E-4, -2.0182056721659587E-5, 7.231106427389245E-5};
		e.expCF(dt, Q.Q, v);
		compare(v, solv, 1e-7);
		
		Q.initA2Q(N);
		double [] A2 = readDoubleMatrix(TEST1_DIR + "A2.txt");
		compare(A2, Q.A2);
		
//		double [] A2Q = readDoubleMatrix(TEST1_DIR + "A2Q.txt");
//		compare(A2Q, Q.A2Q);

		
		COMPLEX [][] w = new COMPLEX[e.ci_real.length][v.length];
		for (int i = 0; i < v.length; i++) {
			for (int j = 0; j < e.ci_real.length; j++) {
				w[j][i] = new COMPLEX();
			}
		}
		
		
		
//	    w = np.zeros([len(c_real),len(v)],dtype=np.complex)
		Arrays.fill(v, 1.0);
		
		v = new double[]{0.006463720994063693, -0.0025598894967633185, -0.011670084506784329, 0.006992547947822397, 0.008307368929750345, -0.009624494344782747, -0.003891770417487691, 0.009999979361246672, -3.164138031256699E-4, -0.008415978246663592, 0.003324991505858777, 0.0057013946826780845, -0.004708080014436136, -0.0028239150287447754, 0.004625163434420616, 5.276322657217052E-4, -0.0036109917030115145, 8.51152557681401E-4, 0.002282182301589475, -0.0013607669299408828, -0.0011063596425523177, 0.0012769476496068285, 3.077497378592821E-4, -9.215487768382256E-4, 1.0234023463673997E-4, 5.388342382900541E-4, -2.3224876169316427E-4, -2.5401044059004865E-4, 2.1803729621996514E-4, 9.22966179980842E-5, -1.679166495819626E-4, -2.0182056721659587E-5, 7.231106427389245E-5};

		double [] LL = Q.Q;
				
		COMPLEX [] vi = new COMPLEX[v.length];
		for (int i = 0; i < vi.length; i++) {
			vi[i] = new COMPLEX(v[i]/dt,0);
		}		
		
		COMPLEX [] AA = new COMPLEX[LL.length];
		// TODO: initialise AA
		for (int i = 0; i < AA.length; i++) {
			AA[i] = new COMPLEX();
		}		
		
		COMPLEX [] A2c = new COMPLEX[N];
		//QMatrix.initA2(N);
		for (int i = 0; i < v.length; i++) {
			double sum = 0;
			for (int k = 0; k < N; k++) {
				sum += v[k] * QMatrix.A2[i*N + k];
			}
			A2c[i] = new COMPLEX(sum/dt, 0);
		}		
		
		COMPLEX [] z_iExp = readComplexArray(TEST3_DIR + "z_i.txt");
		
		for (int i = 0; i < MatrixExponentiator.ci_real.length; i++) {
			System.err.println("Round " + i);
			// vi = np.sum(w,axis=0)
			double zi_real = MatrixExponentiator.zi_real[i]/dt; 
			double zi_imag = MatrixExponentiator.zi_imag[i]/dt;
			
			assertEquals(z_iExp[i].m_fRe, zi_real, 2e-13);
			assertEquals(z_iExp[i].m_fIm, zi_imag, 2e-13);
			
			int x = 0;
			for (int j = 0; j < AA.length; j++) {
				//for (int k = 0; k < N; k++) {
					AA[x].m_fRe = Q.A2Q[x] - zi_real * Q.A2[x];
					AA[x].m_fIm =          - zi_imag * Q.A2[x];
					x++;
				//}
			}
			
			for (int j = 0; j < v.length; j++) {
				double sum = 0;
				for (int k = 0; k < N; k++) {
					sum += v[k] * QMatrix.A2[j*N + k];
				}
				A2c[j] = new COMPLEX(sum/dt, 0);
			}		

			double [] LLExp = readDoubleMatrix(TEST3_DIR + "LL" + i +".txt");
			compare(LLExp, LL);
			COMPLEX [] viExp = readComplexArray(TEST3_DIR + "c_i" + i +".txt");
			compare(viExp, vi);
			COMPLEX [] A2cExp = readComplexArray(TEST3_DIR + "A2c" + i +".txt");
			compare(A2cExp, A2c);
			COMPLEX [] AAExp = readComplexMatrix(TEST3_DIR + "AA" + i +".txt");
			compare(AAExp, AA);
			
			COMPLEX z = new COMPLEX(zi_real, zi_imag);
			
			// e.fasterSolver(LLExp, AAExp, A2cExp, vi, z);
			e.getValuesUBand(LL, AA, A2c, vi, z, e.out1);
			
			COMPLEX [] array = readComplexArray(TEST3_DIR + "uband0" + i +".txt");
			compareShortest(array, e.out1[0]);
			array = readComplexArray(TEST3_DIR + "uband1" + i +".txt");
			compareShortest(array, e.out1[1]);
			array = readComplexArray(TEST3_DIR + "uband2" + i +".txt");
			compareShortest(array, e.out1[2]);
			array = readComplexArray(TEST3_DIR + "uband3" + i +".txt");
			compareShortest(array, e.out1[3]);
			
			e.ubandSolver(e.out1, e.solvEven, N);
			array = readComplexArray(TEST3_DIR + "sol_even" + i +".txt");
			compareShortest(array, e.solvEven);

			e.getValuesTridag(AA, A2c, e.out2);
			array = readComplexArray(TEST3_DIR + "triband0" + i +".txt");
			compareShortest(array, e.out2[0]);
			array = readComplexArray(TEST3_DIR + "triband1" + i +".txt");
			compareShortest(array, e.out2[1]);
			array = readComplexArray(TEST3_DIR + "triband2" + i +".txt");
			compareShortest(array, e.out2[2]);
			array = readComplexArray(TEST3_DIR + "triband3" + i +".txt");
			compareShortest(array, e.out2[3]);
			
			e.triDiagSolver(e.out2, e.solvOdd, N);
			array = readComplexArray(TEST3_DIR + "sol_odd" + i +".txt");
			compareShortest(array, e.solvOdd, 1e-6);

			for (int k = A2c.length - 1; k >= 1; k -= 2) {
				e.solvEven[k].m_fRe = e.solvEven[k/2].m_fRe;
				e.solvEven[k].m_fIm = e.solvEven[k/2].m_fIm;
				e.solvEven[k-1].m_fRe = e.solvOdd[k/2-1].m_fRe;
				e.solvEven[k-1].m_fIm = e.solvOdd[k/2-1].m_fIm;
			}

			COMPLEX [] xExp = readComplexArray(TEST3_DIR + "sol" + i +".txt");
			compare(xExp, e.solvEven, 1e-8);

			
			for (int j = 0; j < N; j++) {
				w[i][j].mul(e.solvEven[j], e.ci_real[i], e.ci_imag[i]);
			}
//			COMPLEX [] wExpected = readComplexMatrix(TEST2_DIR + "w" + i +".txt");
//			compare(wExpected, w, 1e-7);
		}
		
		w = new COMPLEX[v.length][e.ci_real.length];
		for (int i = 0; i < v.length; i++) {
			for (int j = 0; j < e.ci_real.length; j++) {
				w[i][j] = new COMPLEX();
			}
		}
//		e.expCF(v, Q, dt, w);
		//COMPLEX [] wExpected = readComplexMatrix(TEST1_DIR + "w.txt");
		//compare(wExpected, w, 1e-7);

		
		

//		e.expCF(v, Q, dt);
//		compare(v, solv, 1e-7);

	}
	
	private void compare(COMPLEX[] out, COMPLEX[] d) {
		compare(out, d, 1e-10);
	}
	
	private void compare(COMPLEX[] out, COMPLEX[] d, double epsilon) {
		assertEquals(out.length, d.length);
		for (int i = 0; i < out.length; i++) {
			if (Math.abs(d[i].m_fRe - out[i].m_fRe) > epsilon ||
				Math.abs(d[i].m_fIm - out[i].m_fIm) > epsilon) {
				System.err.println(i + " " + d[i] + " != " + out[i]);
			}
			assertEquals(d[i].m_fRe, out[i].m_fRe, epsilon);
			assertEquals(d[i].m_fIm, out[i].m_fIm, epsilon);
		}
	}

	private void compareShortest(COMPLEX[] out, COMPLEX[] d) {
		compareShortest(out, d, 1e-10);
	}
	private void compareShortest(COMPLEX[] out, COMPLEX[] d, double epsilon) {
		int n = Math.min(out.length, d.length);
		for (int i = 0; i < n; i++) {
			if (Math.abs(d[i].m_fRe - out[i].m_fRe) > epsilon ||
				Math.abs(d[i].m_fIm - out[i].m_fIm) > epsilon) {
				System.err.println(i + " " + d[i] + " != " + out[i]);
			}
			assertEquals(d[i].m_fRe, out[i].m_fRe, epsilon);
			assertEquals(d[i].m_fIm, out[i].m_fIm, epsilon);
		}
	}

	private void compare(COMPLEX[] out, COMPLEX[][] d) {
		compare(out, d, 1e-10);
	}
	
	private void compare(COMPLEX[] out, COMPLEX[][] d, double epsilon) {
		int N = d[0].length;
		//assertEquals(out.length, d.length * N);
		for (int i = 0; i < out.length; i++) {
			assertEquals(d[i/N][i%N].m_fRe, out[i].m_fRe, epsilon);
			assertEquals(d[i/N][i%N].m_fIm, out[i].m_fIm, epsilon);
		}
	}

	private void compare(double[] out, double[] d) {
		compare(out, d, 1e-10);
	}
	
	private void compare(double[] out, double[] d, double epsilon) {
		assertEquals(out.length, d.length);
		int n = 0;
		for (int i = 0; i < out.length; i++) {
			if (Math.abs(d[i] - out[i]) > epsilon) {
				System.err.println(i + " " + d[i] + " != " + out[i]);
				n++;
			}
		}
		if (n > 0) {
			System.err.println(n + " mismatches");
		}
		for (int i = 0; i < out.length; i++) {
			assertEquals(d[i], out[i], epsilon);
		}
	}

	private void init(COMPLEX[] array) {
		for (int i = 0; i < array.length; i++) {
			array[i] = new COMPLEX();
		}
	}
	
	
	
}
