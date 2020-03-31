package com.company.hcad_cli;

import java.io.*;
import java.util.*;

import com.opencsv.exceptions.CsvException;
import search.RangeSearcher;
import btree4j.Value;
import btree4j.indexer.BasicIndexQuery;
import btree4j.indexer.BasicIndexQuery.IndexConditionBW;
import btree4j.BTreeException;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.concurrent.ExecutionException;

import com.opencsv.*;
public class Main {
    private static final int refColCount = 67140;
    private static RangeSearcher rangeSearcher;
    private static Scanner sc = new Scanner(System.in);
    static {
        try {
            rangeSearcher = new RangeSearcher("data");
        } catch (BTreeException e) {
            e.printStackTrace();
        }
    }

    private static void insertMatrix() throws IOException, BTreeException {
        int mappingCount = 0;
        File file = new File("gene_name.tsv");
        String strLine;
        BufferedReader bufferedReader = new BufferedReader(new FileReader(file));
        while (null != (strLine = bufferedReader.readLine())) {
            mappingCount++;
        }

        file = new File("count.txt");
        bufferedReader = new BufferedReader(new FileReader(file));
        int curRow = 0;
        while (null != (strLine = bufferedReader.readLine())) {
            curRow = Integer.valueOf(strLine);
        }

        int mapping[] = new int[mappingCount];
        String colName[] = new String[refColCount];
        double data[] = new double[refColCount];
        file = new File("ref.txt");
        bufferedReader = new BufferedReader(new FileReader(file));
        strLine = null;
        HashMap<String, Integer> refTable = new HashMap<String, Integer>();
        long start = System.currentTimeMillis();
        int lineCount = 1;
        while (null != (strLine = bufferedReader.readLine())) {
            refTable.put(strLine, lineCount - 1); //制作参考表
            colName[lineCount - 1] = strLine; //制作列名，注意数组下标从0开始
            if (curRow == 0) rangeSearcher.addColumn(strLine); //如果是空白数据库，则添加column
            lineCount++;
        }

        file = new File("gene_name.tsv");
        bufferedReader = new BufferedReader(new FileReader(file));
        lineCount = 1;
        while (null != (strLine = bufferedReader.readLine())) {
            if (refTable.get(strLine) != null)
                mapping[lineCount - 1] = refTable.get(strLine); //制作映射，这里的第一列对应refTable的第？列
            else mapping[lineCount - 1] = -1;
            lineCount++;
        }
        refTable = null; //完成历史使命了

        file = new File("matrix.mtx");
        bufferedReader = new BufferedReader(new FileReader(file));
        lineCount = 1;
        int lastRow = 0;
        //int cellCount=1;
        int cellCount = curRow + 1; //是这样的，如果没有细胞，则从第一个开始；若已有n个，则从第n+1个开始。
        int col = 0, row = 0;
        double value = 0;
        while (null != (strLine = bufferedReader.readLine())) {
            //strLine为本行数据
            String[] strArr = strLine.split(" ");
            row = Integer.valueOf(strArr[0]); //行，列
            col = Integer.valueOf(strArr[1]);
            col = mapping[col - 1]; //偏移一下
            value = Double.valueOf(strArr[2]);
            if (lineCount != 1 && row != lastRow) {
                //说明上一个细胞的数据已经输入完毕
                //插入
                rangeSearcher.insert(cellCount, colName, data, refColCount); //插一行进去
                if (cellCount % 100 == 0) System.out.println("当前插入第" + cellCount + "个");
                if (cellCount % 1000 == 0) {
                    rangeSearcher.flush();
                }
                cellCount++;
                Arrays.fill(data, 0); //清零数据，以免出问题
                if (col != -1) data[col] = value; //新细胞
            } else {
                if (col != -1) data[col] = value;
            }
            lastRow = row;
            lineCount++;
        }
        rangeSearcher.insert(cellCount, colName, data, refColCount);
        rangeSearcher.flush();
        file = new File("count.txt");
        BufferedWriter out = new BufferedWriter(new FileWriter(file));
        out.write(String.valueOf(cellCount));
        out.close();
    }

    private static void insertMetadata() throws ClassNotFoundException, SQLException, IOException, CsvException {
        Class.forName("org.sqlite.JDBC");// 加载驱动,连接sqlite的jdbc
        Connection connection = DriverManager.getConnection("jdbc:sqlite:metadata.db");//连接数据库metadata.db
        Statement statement = connection.createStatement();   //创建连接对象，是Java的一个操作数据库的重要接口
        String fileName = "metadata.tsv";
        InputStreamReader is = new InputStreamReader(new FileInputStream(fileName), "utf-8");
        CSVParser csvParser = new CSVParserBuilder().withSeparator('\t').build();
        CSVReader reader = new CSVReaderBuilder(is).withCSVParser(csvParser).build();
        List<String[]> strings = reader.readAll();
        int count=0;
        for (String[] strs : strings) {
            if (count == 0) {
                //第一行是字段，不用管
            } else {
                statement.executeUpdate("INSERT INTO metadata (organ,region,subregion,sample_status,cell_ID,donor_ID,donor_gender,donor_age) " +
                        "VALUES ('"+strs[0]+"','"+strs[1]+"','"+strs[2]+"','"+strs[3]+"','"+strs[4]+"','"+strs[5]+"','"+strs[6]+"','"+strs[7]+"')");
            }
            count++;
        }
//        System.out.println("insert completed");
        reader.close();
        is.close();
        connection.close();//关闭数据库连
    }
    private static void queryByMatrix() throws InterruptedException, ExecutionException, BTreeException {
        int conditionCount = 0;
        double a, b;
        System.out.println("请输入条件组数：");
        conditionCount = sc.nextInt();
        sc.nextLine();
        BasicIndexQuery[] conditions = new BasicIndexQuery[conditionCount];
        String[] colNames = new String[conditionCount];
        for (int i = 1; i <= conditionCount; i++) {
            System.out.println("请输入列名：");
            colNames[i - 1] = sc.nextLine();
            System.out.println("请输入条件（下界 上界）：");
            a = sc.nextDouble();
            b = sc.nextDouble();
            sc.nextLine();
            conditions[i - 1] = new IndexConditionBW(new Value(a), new Value(b));
        }
        Set<Long> resultSet = rangeSearcher.rangeSearch(colNames, conditions, conditionCount);
    }
    private static void queryByMetadata() throws ClassNotFoundException, SQLException {
        Class.forName("org.sqlite.JDBC");// 加载驱动,连接sqlite的jdbc
        Connection connection = DriverManager.getConnection("jdbc:sqlite:metadata.db");//连接数据库metadata.db
        Statement statement = connection.createStatement();   //创建连接对象，是Java的一个操作数据库的重要接口
        String query="SELECT * FROM metadata WHERE ";
        String colName=null;
        String condition=null;
        System.out.println("请输入条件组数：");
        int conditionCount = sc.nextInt();
        sc.nextLine();
        for (int i = 1; i <= conditionCount; i++) {
            System.out.println("请输入列名：");
            colName=sc.nextLine();
            System.out.println("请输入条件（只能精确查询）：");
            condition=sc.nextLine();
            if (i==1) query=query+colName+"='"+condition+"' ";
            else query=query+"AND "+colName+"='"+condition+"' ";
        }
        System.out.println(query);
        statement.executeQuery(query);
        connection.close();//关闭数据库连
    }
    public static void main(String[] args) throws BTreeException, IOException, SQLException, ClassNotFoundException, ExecutionException, InterruptedException, CsvException {
        System.out.println("请输入要进行的操作（1：插入（metadata和matrix）；2：按metadata信息查询；3：按matrix信息查询）：");
        int operation = sc.nextInt();
        switch (operation) {
            case 1:
                insertMatrix();
                insertMetadata();
                break;
            case 2:
                queryByMetadata();
                break;
            case 3:
               queryByMatrix();
                break;
            default:
                break;
        }
    }
}
