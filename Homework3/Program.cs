using System;
using System.CodeDom.Compiler;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MultipleSequenceAligment
{

    struct matrixStruct
    {
        public int val; //which score will be put in the box
        public string Up, Left, Diag; //for scores
        public char proteinUP; //for char up
        public char proteinLeft; //for char left
    }

    struct similarityStruct
    {
        public double similarityValue;
        public string seq1_name;
        public string seq2_name;
        public bool isMerged;
        public string align1_name;
        public string align2_name;
    }

    struct GuideTree
    {
        public string seq1;
        public string seq2;
        public double val;
        public string align1_name;
        public string align2_name;
    }

    class Program
    {

        public static List<GuideTree> guideTree = new List<GuideTree>();
        public static string[] temp2 = new string[4];
        static void Main(string[] args)
        {
            Console.Write("Enter gap penalty:");
            int gapPenalty = Convert.ToInt32(Console.ReadLine());
            InitializationStep(gapPenalty);
            MergeAlignment(temp2, gapPenalty);
            Console.ReadKey();

        }

        public static void InitializationStep(int gapPenalty)
        {

            string text = System.IO.File.ReadAllText(@"C:\Users\TOSHIBA\Desktop\input2.txt");
            string[] proteins = text.Trim().Split('>');
            ArrayList arr = new ArrayList();
            List<similarityStruct> similarity = new List<similarityStruct>();
            int counter = 0;


            for (int k = 0; k < proteins.Length; k++)
            {
                string[] proteinSeq1 = proteins[k].Trim().Split('\n');
                if (proteinSeq1.Length > 1)
                {
                    string alignSeq1 = proteinSeq1[1].Trim();
                    for (int l = k + 1; l < proteins.Length; l++)
                    {
                        string[] proteinSeq2 = proteins[l].Trim().Split('\n');

                        if (proteinSeq2.Length > 1)
                        {
                            string alignSeq2 = proteinSeq2[1].Trim();
                            matrixStruct[,] scoringMatrix =
                                new matrixStruct[alignSeq1.Length + 1, alignSeq2.Length + 1];


                            //Initialization Step - filled with gap penalty for the first row and the first column of matrix
                            for (int i = 1; i < scoringMatrix.GetLength(0); i++)
                            {
                                scoringMatrix[i, 0].val = scoringMatrix[i - 1, 0].val + gapPenalty;
                                scoringMatrix[i, 0].Up = "Up";

                                scoringMatrix[i, 0].proteinLeft = Convert.ToChar(alignSeq1.Substring(i - 1, 1));
                                scoringMatrix[i, 0].proteinUP = '-';
                            }

                            for (int j = 1; j < scoringMatrix.GetLength(1); j++)
                            {
                                scoringMatrix[0, j].val = scoringMatrix[0, j - 1].val + gapPenalty;
                                scoringMatrix[0, j].Left = "Left";

                                scoringMatrix[0, j].proteinUP = Convert.ToChar(alignSeq2.Substring(j - 1, 1));
                                scoringMatrix[0, j].proteinLeft = '-';

                            }

                            //Matrix Fill Step
                            for (int i = 1; i <= alignSeq1.Length; i++)
                            {
                                for (int j = 1; j <= alignSeq2.Length; j++)
                                {

                                    char c1 = Convert.ToChar(alignSeq1.Substring(i - 1, 1));
                                    char c2 = Convert.ToChar(alignSeq2.Substring(j - 1, 1));
                                    int score = BlosumMatrix(c1, c2);

                                    int scoreDiag = 0;
                                    if (alignSeq1.Substring(i - 1, 1) == alignSeq2.Substring(j - 1, 1))
                                    {
                                        scoreDiag = scoringMatrix[i - 1, j - 1].val + score;
                                        scoringMatrix[i, j].proteinUP = Convert.ToChar(alignSeq1.Substring(i - 1, 1));
                                        scoringMatrix[i, j].proteinLeft = Convert.ToChar(alignSeq2.Substring(j - 1, 1));
                                    }
                                    else
                                    {
                                        scoreDiag = scoringMatrix[i - 1, j - 1].val + score;
                                        scoringMatrix[i, j].proteinLeft = Convert.ToChar(alignSeq1.Substring(i - 1, 1));
                                        scoringMatrix[i, j].proteinUP = Convert.ToChar(alignSeq2.Substring(j - 1, 1));
                                    }
                                    int scoreLeft = scoringMatrix[i, j - 1].val + gapPenalty;
                                    int scoreUp = scoringMatrix[i - 1, j].val + gapPenalty;
                                    int maxScore = Math.Max(Math.Max(scoreDiag, scoreLeft), scoreUp);
                                    if (maxScore == scoreLeft)
                                    {
                                        scoringMatrix[i, j].val = maxScore;
                                        scoringMatrix[i, j].Left = "Left";
                                    }

                                    if (maxScore == scoreUp)
                                    {
                                        scoringMatrix[i, j].val = maxScore;
                                        scoringMatrix[i, j].Up = "Up";
                                    }

                                    if (maxScore == scoreDiag)
                                    {
                                        scoringMatrix[i, j].val = maxScore;
                                        scoringMatrix[i, j].Diag = "Diag";
                                    }
                                }
                            }
                            //Traceback Step
                            string seq1 = "";
                            string seq2 = "";
                            int smCount = 0;
                            // int exactMatchScore = 0;
                            int m = scoringMatrix.GetLength(0) - 1;
                            int n = scoringMatrix.GetLength(1) - 1;
                            while (m >= 0 && n >= 0)
                            {
                                if (scoringMatrix[m, n].Left == "Left")
                                {
                                    seq2 += scoringMatrix[m, n].proteinUP;
                                    seq1 += '-';
                                    n--;
                                }
                                else if (scoringMatrix[m, n].Up == "Up")
                                {
                                    seq1 += scoringMatrix[m, n].proteinLeft;
                                    seq2 += '-';
                                    m--;

                                }
                                else if (scoringMatrix[m, n].Diag == "Diag")
                                {
                                    seq1 += scoringMatrix[m, n].proteinLeft;
                                    seq2 += scoringMatrix[m, n].proteinUP;
                                    m--;
                                    n--;
                                    if (scoringMatrix[m, n].proteinLeft == scoringMatrix[m, n].proteinUP)
                                    {
                                        smCount++;
                                    }
                                }
                                else
                                {
                                    break;
                                }
                            }
                            similarityStruct sm = new similarityStruct();
                            sm.seq1_name = proteinSeq1[0].Trim().ToString();
                            sm.seq2_name = proteinSeq2[0].Trim().ToString();
                            sm.align1_name = seq1;
                            sm.align2_name = seq2;
                            // sm.similarityValue = (double)scoringMatrix[scoringMatrix.GetLength(0)-1,scoringMatrix.GetLength(1)-1].val/(Double)seq1.Length;//burayı sequence length e bölcez
                            sm.similarityValue = (Double)smCount / (Double)seq1.Length;
                            similarity.Add(sm);


                            if (counter < 4)
                            {
                                if (l == 1)
                                {
                                    temp2[counter] = seq1;
                                    temp2[counter] = seq2;
                                }
                                else
                                {
                                    temp2[counter] = seq1;
                                }

                            }

                            // seq1 += "  <"+proteinSeq1[0].Trim().ToString();
                            //   seq2 +=  "  <"+proteinSeq2[0].Trim().ToString();
                            //print sequences

                            counter++;

                            PrintSeq(seq1, seq2);

                            Console.WriteLine();
                        }
                    }
                }
            }
            //guide tree
            GuideTree(similarity);

        }
        //initialization Blosum62 matrix
        public static int BlosumMatrix(char c1, char c2)
        {
            int val = 0;
            string[] Blosum62 = System.IO.File.ReadAllText(@"C:\Users\TOSHIBA\Desktop\Blosum62.txt").Trim()
                .Split('\n');
            string[] line = Blosum62[0].Trim().Split(' ');

            int size = 0;

            for (int i = 0; i < line.Length; i++)
            {
                if (line[i] != "")
                {
                    size++;
                }
            }

            string[,] blosumArr = new string[size + 1, size + 1];
            int count = 0;
            string lineClearfy = "";
            string[] Letter = null;
            for (int i = 0; i < Blosum62.Length; i++)
            {
                line = Blosum62[i].Trim().Split(' ');

                for (int m = 0; m < line.Length; m++)
                {
                    if (line[m] != "")
                    {
                        lineClearfy += line[m] + ", ";
                    }
                }

                line = null;

                line = lineClearfy.Trim().Split(',');

                for (int k = 0; k < blosumArr.GetLength(1); k++)
                {
                    if (i != 0)
                    {
                        blosumArr[count, k] = line[k + 1];
                    }
                }

                if (i != 0)
                {
                    count++;
                }

                if (i == 0)
                {
                    Letter = lineClearfy.Split(',');
                }

                lineClearfy = "";
            }

            int cx = 0;
            int cy = 0;

            bool flag1 = true;
            bool flag2 = true;

            for (int j = 0; j < blosumArr.GetLength(1); j++)
            {
                if (Letter[j].Trim().ToString() == c1.ToString().Trim() && flag1)
                {
                    flag1 = false;
                    cx = j;
                }

                if (Letter[j].Trim().ToString() == c2.ToString().Trim() && flag2)
                {
                    flag2 = false;
                    cy = j;
                }
            }

            val = Convert.ToInt32(blosumArr[cy, cx]);

            return val;
        }

        public static void PrintSeq(String seq1, string seq2)
        {

            for (int i = seq1.Length - 1; i >= 0; i--)
            {
                Console.Write(seq1[i]);
            }

            Console.WriteLine();
            for (int i = seq2.Length - 1; i >= 0; i--)
            {
                Console.Write(seq2[i]);
            }

            Console.WriteLine();
        }

        public static GuideTree GuideTree(List<similarityStruct> similarity)
        {
            List<String> seq = new List<String>();
            List<similarityStruct> new_sm = new List<similarityStruct>();
            foreach (var item in similarity)
            {
                if (!seq.Contains(item.seq1_name))
                {
                    seq.Add(item.seq1_name);
                    seq.Add(item.align1_name);
                }
                else if (!seq.Contains(item.seq2_name))
                {
                    seq.Add(item.seq2_name);
                    seq.Add(item.align2_name);
                }
            }

            int length = seq.Count();

            similarityStruct sm = new similarityStruct();
            double val = 0;
            foreach (var item in similarity)
            {
                if (item.similarityValue > val)
                {
                    val = item.similarityValue;
                    sm = item;
                }
            }

            new_sm.Add(sm);
            GuideTree gd;
            gd = new GuideTree();
            gd.seq1 = sm.seq1_name;
            gd.seq2 = sm.seq2_name;
            gd.align1_name = sm.align1_name;
            gd.align2_name = sm.align2_name;
            gd.val = sm.similarityValue;
            guideTree.Add(gd);

            List<similarityStruct> new_similarity = new List<similarityStruct>();
            foreach (var item in similarity)
            {
                if (!sm.Equals(item))
                {
                    new_similarity.Add(item);
                }
            }

            bool flag1 = false;
            bool flag2 = false;
            bool flag3 = false;
            int count = 0;
            string s1 = "";
            string s2 = "";
            List<similarityStruct> temp = new List<similarityStruct>();

            for (int i = 0; i < guideTree.Count; i++)
            {
                double x = 0;
                for (int j = 0; j < new_similarity.Count; j++)
                {
                    if (new_similarity.ElementAt(j).seq2_name != guideTree.ElementAt(i).seq1 && new_similarity
                                                                                                 .ElementAt(j)
                                                                                                 .seq2_name !=
                                                                                             guideTree.ElementAt(i).seq2
                                                                                            && new_similarity
                                                                                                 .ElementAt(j)
                                                                                                 .seq1_name !=
                                                                                             guideTree.ElementAt(i)
                                                                                                 .seq2 && new_similarity
                                                                                                 .ElementAt(j)
                                                                                                 .seq1_name !=
                                                                                             guideTree.ElementAt(i).seq1
                                                                                             && !new_similarity
                                                                                                 .ElementAt(j).isMerged)
                    {
                        similarityStruct a = new similarityStruct();
                        a = new_similarity.ElementAt(j);
                        a.isMerged = true;
                        new_similarity.RemoveAt(j);
                        new_similarity.Add(a);
                        s1 = "";
                        s2 = "";
                        x = 0;
                        count = 0;
                        temp = new List<similarityStruct>();
                    }

                    if (j < new_similarity.Count && s1 == "")
                    {
                        if (guideTree.ElementAt(i).seq1 == new_similarity.ElementAt(j).seq1_name)
                        {
                            x += new_similarity.ElementAt(j).similarityValue;
                            flag1 = true;
                            s1 = new_similarity.ElementAt(j).seq2_name;
                        }
                    }

                    if (j < new_similarity.Count && s1 == "")
                    {
                        if (guideTree.ElementAt(i).seq1 == new_similarity.ElementAt(j).seq2_name)
                        {
                            x += new_similarity.ElementAt(j).similarityValue;
                            flag2 = true;
                            s1 = new_similarity.ElementAt(j).seq1_name;
                        }
                    }

                    if ((flag1 == false && flag2 == true) || (flag2 == false && flag1 == true))
                    {
                        temp.Add(new_similarity.ElementAt(j));
                        new_similarity.RemoveAt(j);
                        flag3 = true;
                        j = -1;
                        count++;
                    }
                    flag1 = false;
                    flag2 = false;

                    if (j < new_similarity.Count && flag3 != true && s2 == "")
                    {
                        if (guideTree.ElementAt(i).seq2 == new_similarity.ElementAt(j).seq1_name)
                        {
                            x += new_similarity.ElementAt(j).similarityValue;
                            flag1 = true;
                            s2 = new_similarity.ElementAt(j).seq2_name;
                        }
                    }

                    if (j < new_similarity.Count && flag3 != true && s2 == "")
                    {
                        if (guideTree.ElementAt(i).seq2 == new_similarity.ElementAt(j).seq2_name)
                        {
                            x += new_similarity.ElementAt(j).similarityValue;
                            flag2 = true;
                            s2 = new_similarity.ElementAt(j).seq1_name;
                        }
                    }

                    if ((flag1 == false && flag2 == true) || (flag2 == false && flag1 == true))
                    {
                        count++;
                        temp.Add(new_similarity.ElementAt(j));
                        new_similarity.RemoveAt(j);
                        j = -1;
                    }

                    if (count == 2 && s1 == s2)
                    {
                        x = x / 2;
                        similarityStruct s = new similarityStruct();
                        s.seq1_name = s1;
                        s.seq2_name = guideTree.ElementAt(i).seq1 + "/" + guideTree.ElementAt(i).seq2;
                        s.similarityValue = x;
                        s.isMerged = true;
                        new_similarity.Add(s);
                        s1 = "";
                        s2 = "";
                        x = 0;
                        count = 0;
                        temp = new List<similarityStruct>();
                    }
                    else if (count == 2)
                    {
                        foreach (var item in temp)
                        {
                            new_similarity.Add(item);
                        }

                        temp = new List<similarityStruct>();
                        count = 0;
                        x = 0;
                    }

                    flag3 = false;
                    flag1 = false;
                    flag2 = false;
                }
            }
            if (new_similarity.Count != 0) return GuideTree(new_similarity);
            else
            {
                gd = new GuideTree();

                if (s1 != "")
                {
                    gd.seq1 = s1;
                }
                else if (s2 != "")
                {
                    gd.seq2 = s2;
                }
                guideTree.Add(gd);
                string str = "";
                foreach (var item in guideTree)
                {
                    if (item.seq1 != null && item.seq2 != null)
                    {
                        str += "[";
                        Console.WriteLine("---------------------");
                    }
                    if (item.seq1 != null)
                    {
                        str += item.seq1 + "/";
                        Console.WriteLine(item.seq1);
                    }
                    if (item.seq2 != null)
                    {
                        str += item.seq2 + "/";
                        Console.WriteLine(item.seq2);

                    }
                    if (item.seq1 != null && item.seq2 != null)
                    {
                        str += "]";
                        Console.WriteLine("---------------------");
                    }
                }



                Console.WriteLine();
                Console.WriteLine(str.ToString());




                return gd;
            }
        }


        struct mergeStruct
        {
            public int val;//which score will be put in the box
            public string Up, Left, Diag;//for scores
            public string proteinUP;//for char up
            public string proteinLeft;//for char left
            public string proteinUP2;//for char up
            public string proteinLeft2;//for char left
        }


        public static void MergeAlignment(string[] seq, int gapPenalty)
        {
            string seq1 = seq[0];
            string seq2 = seq[1];
            string seq3 = seq[2];
            string seq4 = seq[3];

            int seqLength = Math.Min(Math.Min(seq1.Length, seq2.Length), Math.Min(seq3.Length, seq4.Length));

            mergeStruct[,] scoringMatrix = new mergeStruct[seqLength, seqLength];

            for (int i = 1; i < scoringMatrix.GetLength(0); i++)
            {
                scoringMatrix[i, 0].val = scoringMatrix[i - 1, 0].val + 2 * gapPenalty;
                scoringMatrix[i, 0].Up = "Up";

                scoringMatrix[i, 0].proteinLeft = seq1.Substring(i - 1, 1);

                scoringMatrix[i, 0].proteinUP = "-";

            }
            for (int j = 1; j < scoringMatrix.GetLength(1); j++)
            {

                scoringMatrix[0, j].val = scoringMatrix[0, j - 1].val + 2 * gapPenalty;
                scoringMatrix[0, j].Left = "Left";

                scoringMatrix[0, j].proteinUP = (seq1.Substring(j - 1, 1));
                scoringMatrix[0, j].proteinLeft = "-";

            }
            //Matrix Fill Step
            for (int i = 1; i < seqLength; i++)
            {
                for (int j = 1; j < seqLength; j++)
                {

                    int diag = 0;
                    int left = 0;
                    int up = 0;
                    int scoreDiag = 0;
                    int scoreUp = 0;
                    int scoreLeft = 0;

                    if (seq1[j - 1] != '-' && seq2[j - 1] != '-' && seq3[j - 1] != '-' && seq4[j - 1] != '-')
                    {
                        diag += BlosumMatrix(seq1[j - 1], seq3[j - 1]);
                        diag += BlosumMatrix(seq1[j - 1], seq4[j - 1]);
                        diag += BlosumMatrix(seq2[j - 1], seq3[j - 1]);
                        diag += BlosumMatrix(seq2[j - 1], seq4[j - 1]);
                        scoreDiag = scoringMatrix[i - 1, j - 1].val + diag;
                        scoringMatrix[i, j].proteinUP = seq1.Substring(i - 1, 1);
                        scoringMatrix[i, j].proteinUP2 = seq2.Substring(i - 1, 1);
                        scoringMatrix[i, j].proteinLeft = seq3.Substring(j - 1, 1);
                        scoringMatrix[i, j].proteinLeft2 = seq4.Substring(j - 1, 1);
                    }
                    else
                    {
                        diag += BlosumMatrix(seq1[j - 1], seq3[j - 1]);
                        diag += BlosumMatrix(seq1[j - 1], seq4[j - 1]);
                        diag += BlosumMatrix(seq2[j - 1], seq3[j - 1]);
                        diag += BlosumMatrix(seq2[j - 1], seq4[j - 1]);
                        scoreDiag = scoringMatrix[i - 1, j - 1].val + diag;
                        scoringMatrix[i, j].proteinLeft = seq1.Substring(i - 1, 1);
                        scoringMatrix[i, j].proteinLeft2 = seq2.Substring(i - 1, 1);
                        scoringMatrix[i, j].proteinUP = seq3.Substring(j - 1, 1);
                        scoringMatrix[i, j].proteinUP2 = seq4.Substring(j - 1, 1);
                    }
                    scoreLeft = scoringMatrix[i, j - 1].val + gapPenalty + left;
                    scoreUp = scoringMatrix[i - 1, j].val + gapPenalty;

                    int maxScore = Math.Max(Math.Max(scoreDiag, scoreLeft), scoreUp);
                    if (maxScore == scoreLeft)
                    {
                        scoringMatrix[i, j].val = maxScore;
                        scoringMatrix[i, j].Left = "Left";
                    }

                    if (maxScore == scoreUp)
                    {
                        scoringMatrix[i, j].val = maxScore;
                        scoringMatrix[i, j].Up = "Up";
                    }

                    if (maxScore == scoreDiag)
                    {
                        scoringMatrix[i, j].val = maxScore;
                        scoringMatrix[i, j].Diag = "Diag";
                    }
                }
            }
            //Traceback Step
            string trace1 = "";
            string trace2 = "";
            string trace3 = "";
            string trace4 = "";

            int m = scoringMatrix.GetLength(0) - 1;
            int n = scoringMatrix.GetLength(1) - 1;
            while (m >= 0 && n >= 0)
            {
                if (scoringMatrix[m, n].Left == "Left")
                {
                    trace3 += scoringMatrix[m, n].proteinUP;
                    trace4 += scoringMatrix[m, n].proteinUP2;
                    trace1 += '-';
                    trace2 += '-';
                    n--;
                }
                else if (scoringMatrix[m, n].Up == "Up")
                {
                    trace1 += scoringMatrix[m, n].proteinLeft;
                    trace2 += scoringMatrix[m, n].proteinLeft2;
                    trace3 += '-';
                    trace4 += '-';
                    m--;
                }
                else if (scoringMatrix[m, n].Diag == "Diag")
                {
                    trace1 += scoringMatrix[m, n].proteinLeft;
                    trace2 += scoringMatrix[m, n].proteinLeft2;
                    trace3 += scoringMatrix[m, n].proteinUP;
                    trace4 += scoringMatrix[m, n].proteinUP2;
                    m--;
                    n--;
                }
                else
                {
                    break;
                }
            }

            //for (int i = 0; i <= seq1.Length; i++)
            //{
            //    for (int j = 0; j <= seq2.Length; j++)
            //    {
            //        Console.Write(scoringMatrix[i, j].val + " ");
            //    }

            //    Console.WriteLine();
            //}
            //Console.WriteLine();
            Console.WriteLine();

            for (int i = trace1.Length - 1; i >= 0; i--)
            {
                Console.Write(trace1[i]);

            }

            Console.WriteLine();
            for (int i = trace2.Length - 1; i >= 0; i--)
            {
                Console.Write(trace2[i]);

            }
            Console.WriteLine();
            for (int i = trace3.Length - 1; i >= 0; i--)
            {
                Console.Write(trace3[i]);

            }
            Console.WriteLine();
            for (int i = trace4.Length - 1; i >= 0; i--)
            {
                Console.Write(trace4[i]);
            }

            Console.ReadLine();

        }

    }

}


