using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using HTLib2;
using HTLib2.Bioinfo;

namespace proj0
{
    partial class Program
    {
        static HashSet<object> lockeds = new HashSet<object>();
        static void Main(string[] args)
        {
            string[] pdbids = { "a", "b", "c", "d" };
            foreach(string pdbid in pdbids)
            {
                string result_path1 = @"E:\sxa5362\sbnma\" + pdbid + "_sbnma.txt";
                string result_path2 = @"E:\sxa5362\sbnma\" + pdbid + "_nma.txt";
                if (HFile.ExistsAll(result_path1, result_path2) == false)
                {
                    var locked = HFile.LockFile(@"E:\sxa5362\sbnma\lock\" + pdbid);
                    if (locked == null)
                        continue;
                    lockeds.Add(locked);
                    // ....
                    HFile.WriteAllText(result_path1, "abc");
                    HFile.WriteAllText(result_path2, "abc");
                }
            }


            if (HFile.Exists(@"C:\Users\sxa5362\Desktop\test\BF_SBNMA.TXT") == false)
            {
                var f1 = HFile.LockFile(@"C:\Users\sxa5362\Desktop\test\1ubq_minimizedlock.txt");
                var f2 = HFile.LockFile(@"C:\Users\sxa5362\Desktop\test\1ubq_minimizedlock.txt");
                var f3 = HFile.LockFile(@"C:\Users\sxa5362\Desktop\test\lock\1a6g.txt");
                var f4 = HFile.LockFile(@"C:\Users\sxa5362\Desktop\test\lock\1a6g.txt");
            }

            {
                Matlab.Register(@"E:\temp\");

                Tinker.Xyz xyz = Tinker.Xyz.FromFile(@"C:\Users\sxa5362\Desktop\test\1ubq_minimized.xyz", false, Tinker.Xyz.Atom.Format.defformat_digit06);
                Tinker.Prm prm = Tinker.Prm.FromFile(@"C:\Users\sxa5362\Desktop\test\charmm22.prm");

                var hess = Tinker.Run.Testhess
                    ( "\""+@"C:\Program Files\Tinker\bin-win64-8.5.3\testhess.exe"+ "\""
                    , xyz
                    , prm
                    , @"E:\temp\"
                    , digits: 20
                    );

                var hessinfo_tinker = Hess.HessInfo.FromTinker(xyz, prm, hess.hess);
                var modes_tinker = hessinfo_tinker.GetModesMassReduced();


                var univ = Universe.BuilderTinker.Build(xyz, prm);
                var sbnma_hessinfo = Hess.GetHessSbNMA(univ, univ.GetCoords(), 12);
            }

            string pdbPath = @"C:\Users\Anthony\Desktop\PdbForStem\";
            Matlab.Register(@"C:\temp\");
            List<string> pdbList = pdbListGen(pdbPath);

            //for (int i = 0; i < pdbList.Count(); i++)
            //{
            //    mainCalc(pdbPath, pdbList[i]);
            //}

            mainCalc(pdbPath, "1AAC.pdb");
        }

        static void mainCalc(string pdbPath, string pdbId)
        {
            Vector[] caCoords = caCoordsGen(pdbPath, pdbId).Item2;
            double[,] cacaDistMatrix = cacaDistMGen(caCoords);
            Matlab.PutMatrix("cacaDistMatrix", cacaDistMatrix);
            double[] anmBfactors = anmBfactorsGen(caCoords, cacaDistMatrix);
            double[] gnmBfactors = gnmBfactorsGen(caCoords, cacaDistMatrix);
            double[] stemBfactors = stemBfactorsGen(caCoords);
            double[] expBfactors = expBfactorsGen(caCoordsGen(pdbPath, pdbId).Item1);
            

            corrCoeff(pdbId, anmBfactors, gnmBfactors, stemBfactors, expBfactors);
        }

        public static List<string> pdbListGen (string pdbPath)
        {
            string line;
            List<string> pdbList = new List<string>();

            try
            {
                StreamReader sr = new StreamReader(pdbPath + "list.txt");
                line = sr.ReadLine();
                while (line != null)
                {
                    pdbList.Add(line);
                    line = sr.ReadLine();
                }
                sr.Close();
            }
            catch (Exception e)
            {
                Console.WriteLine("Exception: " + e.Message);
            }

            return pdbList;
        }

        public static Tuple<Pdb.Atom[], Vector[]> caCoordsGen(string pdbPath, string pdbName)
        {
            Pdb pdb = Pdb.FromFile(pdbPath + pdbName);
            var caAtoms = pdb.atoms.SelectByName("CA").ToArray();
            caAtoms = caAtoms.SelectByDefault().ToArray();
            Vector[] caCoords = caAtoms.ListCoord().ToArray();
            return Tuple.Create(caAtoms, caCoords);
        }

        public static double[,] cacaDistMGen(Vector[] caCoords)
        {
            double[,] cacaDistMatrix = new double[caCoords.Count(), caCoords.Count()];

            for (int i = 0; i < caCoords.Count(); i++)
            {
                cacaDistMatrix[i, i] = 0;
            }
            for (int i = 0; i < (caCoords.Count() - 1); i++)
            {
                for (int j = i + 1; j < caCoords.Count(); j++)
                {
                    double dist2
                        = ((caCoords[i][0] - caCoords[j][0]) * (caCoords[i][0] - caCoords[j][0]))
                        + ((caCoords[i][1] - caCoords[j][1]) * (caCoords[i][1] - caCoords[j][1]))
                        + ((caCoords[i][2] - caCoords[j][2]) * (caCoords[i][2] - caCoords[j][2]));
                    double dist = Math.Sqrt(dist2);
                    cacaDistMatrix[j, i] = cacaDistMatrix[i, j] = dist;
                }
            }

            return cacaDistMatrix;
        }

        public static double[,] anmHMGen(Vector[] caCoords, double[,] cacaDistMatrix)
        {
            double[,] anmHMatrix = new double[caCoords.Count() * 3, caCoords.Count() * 3];
            int hLength = anmHMatrix.GetLength(0);

            for (int i = 0; i < hLength - 1; i++)
            {
                for (int j = i + 1; j < hLength; j++)
                {
                    if (!(((i % 3 == 0) && ((j == i + 1) || (j == i + 2))) || ((i % 3 == 1) && (j == i + 1))))
                    {
                        if (cacaDistMatrix[i / 3, j / 3] > 13)
                        {
                            anmHMatrix[i, j] = 0;
                            anmHMatrix[j, i] = anmHMatrix[i, j];
                        }
                        else
                        {
                            anmHMatrix[i, j] = (-1) * (caCoords[j / 3][i % 3] - caCoords[i / 3][i % 3]) * (caCoords[j / 3][j % 3] - caCoords[i / 3][j % 3]) / (cacaDistMatrix[i / 3, j / 3] * cacaDistMatrix[i / 3, j / 3]);
                            anmHMatrix[j, i] = anmHMatrix[i, j];
                            HDebug.Assert(double.IsNaN(anmHMatrix[i, j]) == false);
                        }
                    }
                }
            }



            for (int i = 0; i < hLength; i++)
            {
                for (int j = i + 1; j < hLength; j++)
                {
                    if ((i % 3 == 0) && (j == i + 1))
                    {
                        for (int k = 0; k < hLength; k++)
                        {
                            if ((k != j) && (k % 3 == 1))
                            {
                                anmHMatrix[i, j] = anmHMatrix[i, j] + (-1) * anmHMatrix[i, k];
                                anmHMatrix[j, i] = anmHMatrix[i, j];
                            }
                        }
                    }
                    else if ((i % 3 == 0) && (j == i + 2))
                    {
                        for (int k = 0; k < hLength; k++)
                        {
                            if ((k != j) && (k % 3 == 2))
                            {
                                anmHMatrix[i, j] = anmHMatrix[i, j] + (-1) * anmHMatrix[i, k];
                                anmHMatrix[j, i] = anmHMatrix[i, j];
                            }
                        }
                    }
                    else if ((i % 3 == 1) && (j == i + 1))
                    {
                        for (int k = 0; k < hLength; k++)
                        {
                            if ((k != j) && (k % 3 == 2))
                            {
                                anmHMatrix[i, j] = anmHMatrix[i, j] + (-1) * anmHMatrix[i, k];
                                anmHMatrix[j, i] = anmHMatrix[i, j];
                            }
                        }
                    }
                }
            }


            for (int i = 0; i < hLength; i++)
            {
                for (int j = 0; j < hLength; j++)
                {
                    if (i == j)
                    {
                        if (j % 3 == 0)
                        {
                            for (int k = 0; k < hLength; k++)
                            {
                                if ((k != j) && (k % 3 == 0))
                                {
                                    anmHMatrix[i, j] = anmHMatrix[i, j] + (-1) * anmHMatrix[i, k];
                                }
                            }
                        }
                        if (j % 3 == 1)
                        {
                            for (int k = 0; k < hLength; k++)
                            {
                                if ((k != j) && (k % 3 == 1))
                                {
                                    anmHMatrix[i, j] = anmHMatrix[i, j] + (-1) * anmHMatrix[i, k];
                                }
                            }
                        }
                        if (j % 3 == 2)
                        {
                            for (int k = 0; k < hLength; k++)
                            {
                                if ((k != j) && (k % 3 == 2))
                                {
                                    anmHMatrix[i, j] = anmHMatrix[i, j] + (-1) * anmHMatrix[i, k];
                                }
                            }
                        }
                    }
                }
            }

            return anmHMatrix;
        }

        public static Tuple<double[], double[,]> EigGen(double[,] anmHMatrix)
        {
            Matlab.PutMatrix("H", anmHMatrix, true);
            Matlab.Execute("[V,D] = eig(H);");
            Matlab.Execute("D = diag(D);");
            double[] D = Matlab.GetVector("D");
            double[,] V = Matlab.GetMatrix("V", true);

            return Tuple.Create(D, V);
        }

        public static double[,] anmInvHMGen(double[,] anmHMatrix)
        {
            double[] D = EigGen(anmHMatrix).Item1;
            double[,] V = EigGen(anmHMatrix).Item2;

            double[,] anmInvHMatrix = new double[anmHMatrix.GetLength(0), anmHMatrix.GetLength(0)];
            int invHLength = anmInvHMatrix.GetLength(0);

            for (int k = 6; k < invHLength; k++)
            {
                for (int i = 0; i < invHLength; i++)
                {
                    for (int j = 0; j < invHLength; j++)
                    {
                        anmInvHMatrix[i, j] = anmInvHMatrix[i, j] + ((V[i, k] / D[k] * V[j, k]));
                    }
                }
            }

            return anmInvHMatrix;
        }

        public static double[] bfactorsGen(double[,] anmInvHMatrix)
        {
            double[] bfactors = new double[anmInvHMatrix.GetLength(0) / 3];

            for (int i = 0; i < bfactors.Length; i++)
            {
                bfactors[i]
                    = (anmInvHMatrix[(3 * i), (3 * i)]
                    + anmInvHMatrix[(3 * i) + 1, (3 * i) + 1]
                    + anmInvHMatrix[(3 * i) + 2, (3 * i) + 2]);
            }

            return bfactors;
        }

        public static double[] anmBfactorsGen(Vector[] caCoords, double[,] cacaDistMatrix)
        {
            double[,] anmHMatrix = anmHMGen(caCoords, cacaDistMatrix);
            double[,] anmInvHMatrix = anmInvHMGen(anmHMatrix);
            double[] anmBfactors = bfactorsGen(anmInvHMatrix);

            return anmBfactors;
        }

        public static double[] gnmBfactorsGen(Vector[] caCoords, double[,] cacaDistMatrix)
        {
            double[,] gnmHMatrix = new double[caCoords.Count(), caCoords.Count()];
            int mLength = gnmHMatrix.GetLength(0);
            double[] gnmBfactors = new double[mLength];
            
            for (int i = 0; i < mLength - 1; i++)
            {
                for (int j = i + 1; j < mLength; j++)
                {
                    if (cacaDistMatrix[i, j] > 13)
                    {
                        gnmHMatrix[i, j] = 0;
                        gnmHMatrix[j, i] = gnmHMatrix[i, j];
                    }
                    else
                    {
                        gnmHMatrix[i, j] = -1;
                        gnmHMatrix[j, i] = gnmHMatrix[i, j];
                    }
                }
            }

            for (int i = 0; i < mLength; i++)
            {
                for (int j = 0; j < mLength; j++)
                {
                    if (i != j)
                    {
                        gnmHMatrix[i, i] = gnmHMatrix[i, i] + gnmHMatrix[i, j];
                    }
                }
                gnmHMatrix[i, i] = (-1) * gnmHMatrix[i, i];
            }
            
            Matlab.PutMatrix("gnmHMatrix", gnmHMatrix, true);
            Matlab.Execute("gnmInvHMatrix = pinv(gnmHMatrix);");
            Matlab.Execute("gnmBfactors = diag(gnmInvHMatrix);");
            gnmBfactors = Matlab.GetVector("gnmBfactors");

            return gnmBfactors;
        }

        public static double[] stemBfactorsGen(Vector[] caCoords)
        {
            Matlab.PutMatrix("caCoords", caCoords);
            Matlab.Execute("[hv1, hv2, hv3, hv4] = STeM(caCoords);");
            Matlab.Execute("stemMatrix = hv1+hv2+hv3+hv4;");
            Matlab.Execute("stemInvHMatrix = pinv(stemMatrix);");
            double[,] stemInvHMatrix = Matlab.GetMatrix("stemInvHMatrix", true);
            double[] stemBfactors = bfactorsGen(stemInvHMatrix);

            return stemBfactors;
        }

        public static double[] expBfactorsGen(Pdb.Atom[] caAtoms)
        {
            double[] expBfactors = caAtoms.ListTempFactor().ToArray();

            return expBfactors;
        }

        public static void corrCoeff(string pdbId, double[] anmBfactors, double[] gnmBfactors, double[] stemBfactors, double[] expBfactors) {

            Matlab.PutVector("anmBfactors", anmBfactors);
            Matlab.PutVector("gnmBfactors", gnmBfactors);
            Matlab.PutVector("stemBfactors", stemBfactors);
            Matlab.PutVector("expBfactors", expBfactors);


            Matlab.Execute("anmCorr = corr(anmBfactors, expBfactors);");
            Matlab.Execute("gnmCorr = corr(gnmBfactors, expBfactors);");
            Matlab.Execute("stemCorr = corr(stemBfactors, expBfactors);");


            double anmCorr = Matlab.GetValue("anmCorr");
            double gnmCorr = Matlab.GetValue("gnmCorr");
            double stemCorr = Matlab.GetValue("stemCorr");

            Console.Write(pdbId + "    ");
            Console.Write( "{0:0.0000} "  , anmCorr );
            Console.Write( "  {0:0.0000} ", gnmCorr );
            Console.Write( "  {0:0.0000}" , stemCorr);
            Console.WriteLine();
        }

        public static void C0V7FileGen(Vector[] caCoords, double[,] V, Pdb.Atom[] caAtoms)
        {
            int size = caCoords.Count();

            for (int i = 0; i < 100; i++)
            {
                Vector[] C0V7 = new Vector[size];
                for (int j = 0; j < size; j++)
                {
                    double[] addVector = new double[3];
                    for (int k = 0; k < 3; k++)
                    {
                        addVector[k] = caCoords[j][k] + V[(3 * j) + k, 7] * i;
                    }
                    C0V7[j] = addVector;
                }
                var caAtoms_C0V7 = caAtoms.CloneByUpdateCoord(C0V7);
                var pdb_C0V7 = Pdb.FromAtoms(caAtoms_C0V7);
                string filename = @"C:\temp\Pymol\1aky_C0V7_" + i + ".pdb";
                pdb_C0V7.ToFile(filename);
            }

            Pdb.Atom[] caAtoms = caAtomsVecGen();
            {
                Vector[] C0 = caCoordsVecGen(caAtoms);

                Vector[] V7 = new Vector[C0.Length];
                for (int i = 0; i < V7.Length; i++)
                    V7[i] = new double[] { 0.5, 1, 1.5 };

                for (int k = 1; k < 10; k++)
                {
                    Vector[] C0V7 = new Vector[C0.Length];
                    for (int i = 0; i < V7.Length; i++)
                        C0V7[i] = C0[i] + V7[i] * k;

                    var caAtoms_C0V7 = caAtoms.CloneByUpdateCoord(C0V7);
                    var pdb_C0V7 = Pdb.FromAtoms(caAtoms_C0V7);
                    string filename = pdbPath;
                    pdb_C0V7.ToFile(filename);
                }
            }
        }
    }
}

namespace proj0
{
    partial class Program
    {
        static long Fib(int init, int n)
        {
            return init + xFib(n);
        }
        static long xFib(int n)
        {
            if (n == 0)
                return 0;
            if (n == 1)
                return 1;
            return xFib(n - 1) + xFib(n - 2);
        }

        static HashSet<object> lockeds = new HashSet<object>();
        static void Main(string[] args)
        {
            string pathbase = @"E:\sxa5362\temp\htna\";
            if (HDirectory.Exists(pathbase) == false)
                HDirectory.CreateDirectory(pathbase);

            for (int init = 0; init < 10; init++)
            {
                var locked = HFile.LockFile(pathbase + init + ".lock");
                if (locked == null)
                {
                    System.Console.WriteLine("Fib(" + init + ") is being handled by other program.");
                    continue;
                }
                lockeds.Add(locked);

                string filepath = pathbase + init + ".log";
                if (HFile.Exists(filepath) == false)
                {
                    List<string> lines = new List<string>();
                    for (int n = 0; n < 40; n++)
                    {
                        long fibn = Fib(init, n);
                        lines.Add("fib(" + n + ") => " + fibn);
                    }
                    HFile.WriteAllLines(filepath, lines);
                    System.Console.WriteLine("Fib(" + init + ") is done.");
                }
                else
                {
                    System.Console.WriteLine("Fib(" + init + ") is skipped.");
                }
            }
        }
    }
}
