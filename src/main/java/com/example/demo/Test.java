package com.example.demo;

import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.databind.DeserializationFeature;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.type.MapType;
import com.fasterxml.jackson.databind.type.TypeFactory;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.OptionalLong;
import java.util.Set;
import java.util.Stack;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.LongStream;

public class Test {
    
    private static Map<String, Node> readJson(){
        ObjectMapper objectMapper = new ObjectMapper();
        objectMapper.configure(DeserializationFeature.FAIL_ON_UNKNOWN_PROPERTIES, false);
        TypeFactory typeFactory = objectMapper.getTypeFactory();
        MapType mapType = typeFactory.constructMapType(Map.class, String.class, Node1.class);
        try {
            return objectMapper.readValue("\n" +
                    "{\n" +
                    "\t\"name\" : {\n" +
                    "\t\t\t\"value1\" : 5,\n" +
                    "\t\t\t\"value\" : 5\n" +
                    "\n" +
                    "\t}\n" +
                    "\n" +
                    "\n" +
                    "}", mapType);
        } catch (JsonProcessingException e) {
            e.printStackTrace();
        }
        
     return new HashMap<>();
    
    }
    private static class Tree{
        long node;
        Tree left;
        Tree right;
    
        public Tree(long node) {
            node = node;
        }
    }
    public static String leftOrRight(long[] arr) {
        Tree root = new Tree(arr[0]);
        int i = 1;
        List<Tree> lastLayers = new ArrayList<>();
        lastLayers.add(root);
        while (i < arr.length){
            List<Tree> tempLayer = new ArrayList<>();
            for (Tree node : lastLayers){
                if(arr[i] != -1){
                    node.left = new Tree(arr[i]);
                    tempLayer.add(node.left);
                }
                i++;
                if(arr[i] != -1){
                    node.right = new Tree(arr[i]);
                    tempLayer.add(node.right);
                }
                i++;
            }
            lastLayers = tempLayer;
        }
        
        long sumLeft = recursiveSum(root.left);
        long sumRight = recursiveSum(root.right);
        return sumLeft == sumRight ? "" : sumLeft < sumRight ? "Left" : "Right";    
    }
    
    private static long recursiveSum(Tree left) {
        if(left == null)
            return 0;
        return left.node + recursiveSum(left.left) + recursiveSum(left.right);
    }
    
    public static int[] twoSum(int[] numbers, int target) {
        int output[] = new int[2];
        for(int i = 0; i < numbers.length; i++){
            output[0] = i + 1;
            int num = target - numbers[i];
            int binarySearched = binarySearch(num, numbers, i - 1, numbers.length);
            if(binarySearched > 0){
                output[1] = binarySearched + 1;
                return output;
            }
        }
        return output;
    }
    
    public static int binarySearch(int num, int[] numbers, int l , int r){
        while(r > l + 1){
            int middle = (l + r) / 2;
            if(numbers[middle] <= num){
                l = middle;
            }
            else{
                r = middle;
            }
        }
        if(numbers[l] == num)
            return l;
        else 
            return -1;
    }
    
    private static class ListNode {
      int val;
      ListNode next;
      ListNode() {}
      ListNode(int val) { val = val; }
      ListNode(int val, ListNode next) { val = val; next = next; }
    }
    
    public static ListNode mergeTwoLists(ListNode l1, ListNode l2) {
        ListNode output = new ListNode(Integer.MIN_VALUE);
        ListNode currNode = output;
        while(l1 != null || l2 != null){
            
            Integer value1 = null;
            if(l1 != null){
                value1 = l1.val;
            }
            Integer value2 = null;
            if(l2 != null){
                value2 = l2.val;
            }
            
            if(value1 == null){
                currNode = addElement(currNode, new ListNode(value2));
                l2 = l2.next;
            }
            else if(value2 == null){
                currNode = addElement(currNode, new ListNode(value1));
                l1 = l1.next;
            }
            else if(value1 <= value2){
                currNode = addElement(currNode, new ListNode(value1));
                l1 = l1.next;
            }
            else {
                currNode = addElement(currNode, new ListNode(value2));
                l2 = l2.next;
            }
        }
        return output;
        
    }
    
    public static ListNode addElement(ListNode output, ListNode elem){
        if(output.val == Integer.MIN_VALUE) {
            output.val = elem.val;
        }
        else {
            output.next = elem;
            output = elem;
        }
        return output;
    }
    
    public static Set<List<Integer>> threeSum(int[] nums) {
        Set<List<Integer>> list = new HashSet<List<Integer>>();
        int[][] table = new int[nums.length][nums.length];
        for(int i = 0; i < nums.length; i++){
            for(int j = 0; j < i; j++){
                table[i][j] = nums[i] + nums[j];
            }
        }
        
        for(int i = 0; i < nums.length; i++){
            
            for(int j = 0 ; j < table.length && j != i ; j++){
                for(int k = 0; k < j && k != i; k++){
                    if(table[j][k] + nums[i] == 0){
                        List<Integer> innerList = new ArrayList<>();
                        innerList.add(nums[i]);
                        innerList.add(nums[j]);
                        innerList.add(nums[k]);
                        list.add(innerList);
                    }
                }
            }
            
            
            
        }
        return list;
        
        
        
    }
    public static int romanToInt(String s) {
        int sum = 0;
        if(s.length() == 0)
            return 0;
        int index = 0;
        Map<String, Integer> map = new HashMap();
        map.put("I", 1);
        map.put("V", 5);
        map.put("X", 10);
        map.put("L", 50);
        map.put("C", 100);
        map.put("D", 500);
        map.put("M", 1000);
        map.put("", 0);
        while(index < s.length()){
            int current = map.get(String.valueOf(s.charAt(index)));
            String nextChar = index + 1 < s.length() ? String.valueOf(s.charAt(index+1)) : "";
            if(map.get(nextChar) > current){
                sum += (map.get(nextChar) - current);
                
            }
            else
            {
                sum += current;
            }
            
            index++;
            
        }
        return sum;
    }
    
    public static int myAtoi(String str) {
        String numberString = str.trim();
        if(numberString.length() == 0)
            return 0;
        int startIndex = 0;
        StringBuilder builder = new StringBuilder();
        while(startIndex <  numberString.length()){
            if(Character.isDigit(numberString.charAt(startIndex)) ||
                    numberString.charAt(startIndex) == '-')
                builder.append(numberString.charAt(startIndex));
            startIndex++;
        }
        numberString = builder.toString();
        startIndex = 0;
        char firstChar = numberString.charAt(0);
        if(firstChar == '-'){
            startIndex++;
        }
        int multiplier = (int)Math.pow(10, numberString.length() - startIndex - 1);
        int sum = 0;
        while(startIndex <  numberString.length()){
            char candidate = numberString.charAt(startIndex++);
            int intChar = Character.getNumericValue(candidate) * multiplier;
            sum += intChar;
            multiplier = multiplier / 10;
        }
        return firstChar == '-' ? sum * -1 : sum;
    }
    
    private static Map<String, String> mapping = new HashMap<>();
    static {
        mapping.put("2", "abc");mapping.put("3", "def");mapping.put("4", "ghi");
        mapping.put("5", "jkl");mapping.put("6", "mno");mapping.put("7", "pqrs");
        mapping.put("8", "tuv");mapping.put("9", "wxyz");
    }
    
    public static String longestPalindrome(String s) {
        int longest =  lps(s, 0 , s.length()-1);
        return "";
    }
    
    //apqrsa
    
    private static int lps(String s, int i, int n) {
        //If single character
        if( i == n)
            return 0;
        //If 2 characters and both are same
        if(i + 1 == n && s.charAt(i) == s.charAt(n)){
            return 2;
        }
        
        //If more than 2 characters and first and last are same
        if(s.charAt(i) == s.charAt(n)){
            return lps(s, i+1, n-1) + 2;
        }
        
        //If first and last characters are not same
        return Math.max(lps(s, i, n-1), lps(s, i+1, n));
    }
    
    public static boolean hasPath(long[][] maze, long n) {
        int mazeSize = (int)n;
        int[][] table = new int[mazeSize][mazeSize];
        //first is blocker
        if(maze[0][0] == 1)
            return false;
        
        //Base, first to first is possible
        table[0][0] = 1;
        
        for (int i = 0; i < maze.length; i++) {
            for (int j = 0; j < maze[0].length; j++) {
                if(maze[i][j] == 1 || (i == 0 && j ==0))
                    continue;
                table[i][j] = getCellValue(table, i-1, j) 
                                + getCellValue(table, i, j - 1);
            }
        }
        
        
        return table[mazeSize - 1][mazeSize - 1] > 0;
    }
    
    
    private static Set<String> validExpressions = new HashSet<String>();
    
    private static void reset() {
        validExpressions.clear();
    }
    
    private static void recurse(
            String s,
            int index,
            int leftRem,
            int rightRem,
            int leftCount,
            int rightCount,
            StringBuilder expression) {
        
        // If we have reached the end of string.
        if (index == s.length()) {
            
            // If the current expression is valid.
            if (leftRem == rightRem) {
                String possibleAnswer = expression.toString();
                validExpressions.add(possibleAnswer);
            }
        } else {
    
            char currentCharacter = s.charAt(index);
            int length = expression.length();
    
            if(currentCharacter == '(' && leftRem > 0){
                recurse(s, index + 1, leftRem - 1, rightRem, leftCount , rightCount, expression);
            }
            else if(currentCharacter == ')' && rightRem > 0){
                recurse(s, index + 1, leftRem, rightRem - 1, leftCount, rightCount, expression);
            }
    
           
            if (currentCharacter == '(') {
                expression.append(currentCharacter);
                recurse(s, index + 1, leftRem, rightRem, leftCount + 1, rightCount, expression);
                expression.deleteCharAt(length);
            }
            else if(currentCharacter == ')') {
                expression.append(currentCharacter);
                if(rightCount < leftCount)
                {
                    recurse(s, index + 1, leftRem, rightRem, leftCount, rightCount + 1, expression);
                }
                expression.deleteCharAt(length);
            }
            else
            {
                // If the current character is neither an opening bracket nor a closing one,
                // simply recurse further by adding it to the expression StringBuilder
                expression.append(currentCharacter);
                recurse(s, index + 1, leftRem, rightRem, leftCount, rightCount, expression);
                expression.deleteCharAt(length);
            }
        }
    }
    
    private static int maxSoFar = 0;
    
    public static int largestSumOfNonAdjacsentNumbers(int[] arr){
        if(arr.length == 0){
            return 0;
        }
        if(arr.length == 2){
            return Math.max(arr[0], arr[1]);
        }
        
        adjacentSum(0, 0, arr);
        return maxSoFar;
    }
    
    public static int search(int[] nums, int target) {
        int min = Math.min(nums[0], nums[nums.length - 1]);
        int max = Math.max(nums[0], nums[nums.length - 1]);
        
        int low = 0;
        int high = nums.length - 1;
        int pivot = 0;
        while(low < high){
            int middle = (low + high)/2;
            if(middle + 1 < nums.length && nums[middle] > nums[middle + 1])
                pivot = middle + 1;
            else if(middle -1 >= 0 && nums[middle] < nums[middle -1])
                pivot = middle;
            if(nums[middle] > nums[low])
                low = middle + 1;
            else
                high = middle;
        }
        
        //search in subarray1
        low = -1;
        high = pivot;
        while(low + 1 < high){
            int middle = (low + high)/2;
            if(nums[middle] == target){
                return middle;
            }
            if(nums[middle] <= target){
                low = middle;
            }
            else{
                high = middle;
            }
        }
        
        low = pivot - 1;
        high = nums.length;
        while(low + 1 < high){
            int middle = (low + high)/2;
            if(nums[middle] == target){
                return middle;
            }
            if(nums[middle] < target){
                low = middle;
            }
            else{
                high = middle;
            }
        }
        
        return -1;
        
    }
    
    public static int divide(int dividend, int divisor) {
        
        if(dividend == 0)
            return 0;
        
        int curr =  0;
        int currBig =  0;
        int count = 0;
        
        int absDividend = Math.abs(dividend);
        int absDivisor = Math.abs(divisor);
        if(absDivisor == 1){
            count = absDividend;
        }
        else if(absDivisor > absDividend){
            count = 0;
        }
        else{
            while(currBig < absDividend){
                curr += absDivisor;
                currBig += absDivisor + absDivisor;
                count = count+2;
            }
            
            while(currBig > absDividend){
                currBig -= absDivisor;
                count--;
            }
        }
        
        
        
        if((dividend < 0 && divisor < 0) || (dividend > 0 && divisor > 0)){
            return count;
        }
        else if(dividend < 0 || divisor < 0){
            return -count;
        }
        
        return count;
        
    }
    private static void adjacentSum(int index, int sum, int[] arr) {
        if(index >= arr.length){
            if(sum > maxSoFar)
                maxSoFar = sum;
            return;
        }
        int curr = arr[index];
        //skip the current one
        adjacentSum(index + 1, sum, arr);
        if(curr > 0){
            sum += curr;
            adjacentSum(index + 2, sum, arr);
            sum -= curr;
        }
    }
    
    public static List<String> removeInvalidParentheses(String s) {
        
        reset();
        int left = 0;
        int right = 0;
        for(int i = 0; i < s.toCharArray().length; i++){
            char ch = s.toCharArray()[i];
            
            if(ch == '('){
                left++;
            }
            else if(ch == ')'){
                if(left > 0)
                    left--;
                else
                    right++;
            }
        }
        
        //left denotes number of left to be removed
        //right denotes number of right to be removed
        recurse(s, 0, left, right, 0, 0, new StringBuilder());
        return new ArrayList(validExpressions);
    }
    
    private static int getCellValue(int[][] table, int x , int y){
        if(x < 0 || y < 0)
            return 0;
        return table[x][y];
    }
    
    
    private static String longestPalinString(String s){
        
        int[][] table = new int[s.length()][s.length()];
        int maxLength = 0;
        int start = 0;
        //All one character strings are palindrom
        for (int i = 0; i < table.length; i++) {
            table[i][i] = 1;
            start = i;
            maxLength = 1;
        }
        
        //2 character strings
        for (int i = 0; i < table.length - 1; i++) {
            if(s.charAt(i) == s.charAt(i+1)){
                table[i][i+1] = 1;
                start = i;
                maxLength = 2;
            }
        }
        
        //More than 2 characters string
        for (int k = 3; k <= table.length; k++) {
            for (int i = 0; i < table.length - k + 1; i++) {
                if(s.charAt(i) == s.charAt(i+k-1) && table[i + 1][i+k-2] == 1){
                    table[i][i+k-1] = 1;
                    start = i;
                    maxLength = k;
                }
            }
        }
        
        return s.substring(start, start+maxLength-1);
    }
    
    
    public static long solution(long[] numbers) {
        OptionalLong optionalLong = LongStream.of(numbers).max();
        return optionalLong.orElse(0l);
        
    }
    
    public static double findMedianSortedArrays(int[] num1, int[] num2) {
        List<Integer> lastVisited = new ArrayList(2);
        int p1 = 0;
        int p2 = 0;
        int totalLength = num1.length + num2.length;
        int visited = 0;
        while(visited <= totalLength / 2){
            if(p1 >= num1.length){
                addToLastVisited(num2[p2], lastVisited);
                p2++;
                visited++;
                continue;
            }
            if(p2 >= num2.length){
                addToLastVisited(num1[p1], lastVisited);
                p1++;
                visited++;
                continue;
            }
            
            if(num1[p1] <= num2[p2]){
                addToLastVisited(num1[p1], lastVisited);
                p1++;
                visited++;
            }
            else{
                addToLastVisited(num2[p2], lastVisited);
                p2++;
                visited++;
            }
            
            
        }
        
        return totalLength % 2 == 0 ? (lastVisited.get(0) + lastVisited.get(1))/2 : lastVisited.size() == 1 ? lastVisited.get(0) : lastVisited.get(1);
    }
    
    private static void addToLastVisited(int val, List<Integer>  lastVisited ){
        if(lastVisited.size() == 2){
            lastVisited.remove(0);
            lastVisited.add(val);
        }
        else{
            lastVisited.add(val);
        }
        
    }
    
    public int lengthOfLongestSubstring(String s) {
        int longestTillNow = 0;
        int startPointer = 0;
        int pointer = 0;
        Set<String> visited = new HashSet<>();
        while(pointer < s.length()){
            if(visited.contains(String.valueOf(s.charAt(pointer)))){
                if(pointer - startPointer > longestTillNow)
                    longestTillNow = pointer - startPointer;
                startPointer += 1;
                pointer = startPointer;
                visited.clear();
            }
            else{
                visited.add(String.valueOf(s.charAt(pointer)));
                pointer++;
            }
        }
        return Math.max(pointer - startPointer, longestTillNow);
    }
    
    
    public static int uniquePathsWithObstacles(int[][] obstacleGrid) {
        int a[][] = new int[obstacleGrid.length][obstacleGrid[0].length];
        if (obstacleGrid[0][0] == 1)
            return 0;
        a[0][0] = 1;
        for (int i = 0; i < obstacleGrid.length; i++) {
            for (int j = 0; j < obstacleGrid[i].length; j++) {
                if(i == 0 && j == 0 )
                    continue;
                if(obstacleGrid[i][j] == 1){
                    a[i][j] = 0;
                }
                else{
                    int leftQualify = i - 1 >= 0 ? a[i-1][j] : 0;
                    int rightQualify = j - 1 >=0 ? a[i][j-1] : 0;
                    a[i][j] = leftQualify + rightQualify;
                }
            }
        }
        return a[obstacleGrid.length - 1][obstacleGrid[0].length - 1];
    }
    
    static long magicLand(int[] Amount, int M, int N) {
        Arrays.sort(Amount);
        long sum = 0;
    
        for (int i = Amount.length - 1; i >= M - 1; i--) {
            sum += (Amount[i] * ( fact(i) / (fact(M - 1) * fact(i - M + 1) )));
        }
        
        return (long)(sum % (Math.pow(10, 9) + 7));
    }
    
    public static List<String> letterCombinations(String digits) {
        List<String> combinations = new ArrayList<>();
        generateCombination("", digits, combinations);
        return combinations;
    }
    
    public static void generateCombination(String currentString, String remainingNumbers, List<String> generatedCombinations){
        if(remainingNumbers.length() == 1){
            String allChars = mapping.get(remainingNumbers);
            for (int i = 0; i < allChars.length(); i++) {
                generatedCombinations.add(currentString + allChars.charAt(i));
            }
            return;
        }
        String allChars = mapping.get(String.valueOf(remainingNumbers.charAt(0)));
        for (int i = 0; i < allChars.length(); i++) {
            generateCombination(currentString + allChars.charAt(i), remainingNumbers.substring(1), generatedCombinations);
        }
    }
    
    static long fact(int num){
        long fact = 1;
        while (num > 1){
            fact *= num--;
        }
        return fact;
    }
    
    public static List<List<Integer>> permute(int[] nums) {
        List<List<Integer>> set = new ArrayList<>();
        
        for(int i = 0; i < nums.length; i++){
            
            List<Integer> remainingNums = new ArrayList();
            for(int j = 0; j < nums.length; j++){
                if(j == i)
                    continue;
                remainingNums.add(nums[j]);
            }
            
            List<Integer> list = new ArrayList();
            list.add(nums[i]);
            recursivePermute(remainingNums , set, list);
        }
        return set;
    }
    
    public static long[] minimumOccurrences(long[] numbers) {
       Map<Long, Integer> numberToCounts = new HashMap<>();
       LongStream.of(numbers).forEach(num -> {
           if(!numberToCounts.containsKey(num))
               numberToCounts.put(num, 0);
           numberToCounts.put(num, numberToCounts.get(num) + 1);
       });
    
       TreeMap<Integer, TreeSet<Long>> countToNumbers = new TreeMap<>();
       for (Map.Entry<Long, Integer> entry : numberToCounts.entrySet()){
           int value = entry.getValue();
           if(countToNumbers.containsKey(value)){
               countToNumbers.get(value).add(entry.getKey());
           }
           else{
               TreeSet<Long> sameOccurrenceNumbers = new TreeSet<>();
               sameOccurrenceNumbers.add(entry.getKey());
               countToNumbers.put(value, sameOccurrenceNumbers);
           }
       }
       if(countToNumbers.isEmpty())
           return new long[0];
       long[] output = new long[countToNumbers.firstEntry().getValue().size() * countToNumbers.firstKey()];
       int index = 0;
       for (Long resultValue : countToNumbers.firstEntry().getValue()){
           for (int i = 0; i < countToNumbers.firstKey(); i++) {
               output[index++] = resultValue;
           }
       }
       return output;
    }
    
    public static String leftOrRightHired(long[] arr) {
        if(arr.length == 0)
            return "";
        Node root = buildTree(arr);
        
        long leftSum = sum(root.left);
        long rightSum = sum(root.right);
        
        return leftSum > rightSum ? "Left" : leftSum == rightSum ? "" : "Right";
    }
    
    private static long sum(Node node){
        if(node == null)
            return 0;
        return node.value + sum(node.left) + sum(node.right);
    }
    
    
    private static class Node1 implements Serializable {
        long value;
        
        Node1(){
            
        }
        
        Node1(long value){
            value = value;
        }
        
    }
    
    private static class Node implements Serializable {
        long value;
        Node left;
        Node right;
        
        Node(long value){
            value = value;
        }
        
    }
    private static Node buildTree(long[] arr) {
        Node root = new Node(arr[0]);
        List<Node> lastLayer = new ArrayList();
        lastLayer.add(root);
        int index = 1;
        while(index < arr.length){
            List<Node> temp = new ArrayList();
            for(Node currNode : lastLayer){
                int leftIndex = index++;
                int rightIndex = index++;
                if(leftIndex < arr.length && arr[leftIndex] != -1){
                    currNode.left = new Node(arr[leftIndex]);
                    temp.add(currNode.left);
                }
                if(rightIndex < arr.length && arr[rightIndex] != -1){
                    currNode.right = new Node(arr[rightIndex]);
                    temp.add(currNode.right);
                }
            }
            lastLayer = temp;
        }
        
        return root;
        
        
        
        
        
        
        
        
    }
    
    public static long maxProfit(long[] prices) {
        if(prices.length < 2)
            return 0;
        long profitSoFar = 0;
        TreeSet<Long> alreadyVisitedBuyingPrices = new TreeSet<>();
        for(int i =0; i < prices.length - 1; i++){
            //If previous stock was in lesser value, we must have got highest profit
            // We need to skip all buying prices higher than that
            if(alreadyVisitedBuyingPrices.size() > 0 &&
                    alreadyVisitedBuyingPrices.last() < prices[i]){
                continue;
            }
            for(int j = i + 1; j < prices.length; j++){
                if(prices[j] - prices[i] > profitSoFar){
                    profitSoFar =  prices[j] - prices[i];
                }
                
            }
            alreadyVisitedBuyingPrices.add(prices[i]);
        }
        return profitSoFar;
        
    }
    
    public List<List<Integer>> permuteUnique(int[] nums) {
        
        List<List<Integer>> finalSet = new ArrayList<List<Integer>>();
        if(nums.length == 1)
        {
            finalSet.add(Arrays.asList(nums[0]));
            return finalSet;
        }
        Map<Integer, Integer> counts = new HashMap();
        
        
        for(int i = 0 ; i < nums.length; i++){
            if(!counts.containsKey(nums[i])){
                counts.put(nums[i], 0);
            }
            counts.put(nums[i], counts.get(nums[i]));
        }
        
        for(java.util.Map.Entry<Integer, Integer> entry : counts.entrySet()){
            int key = entry.getKey();
            List<Integer> currentSet = new ArrayList();
            currentSet.add(key);
            int value = entry.getValue();
            List<Integer> remainingChars =
                    IntStream.of(nums).boxed().collect(Collectors.toList());
            remainingChars.remove(entry.getKey());
            
            permute(finalSet, remainingChars, currentSet);
        }
        
        return finalSet;
        
    }
    
    private void permute(List<List<Integer>> finalSet, List<Integer> remainingChars, List<Integer> currentSet){
        if(remainingChars.size() == 1){
            currentSet.add(remainingChars.get(0));
            if(!finalSet.contains(currentSet)){
                finalSet.add(currentSet);
            }
            return;
        }
        
        for(int i = 0 ; i < remainingChars.size(); i++){
            List<Integer> newCurrentSet = new ArrayList();
            newCurrentSet.addAll(currentSet);
            newCurrentSet.add(remainingChars.get(i));
            
            List<Integer> newRemainingChars = new ArrayList(remainingChars);
            newRemainingChars.remove(i);
            
            permute(finalSet, newRemainingChars, newCurrentSet);
            
        }
        
        
        
    }
    private static void recursivePermute(List<Integer> remainingNumbers, List<List<Integer>> superSet, List<Integer> set){
        if(remainingNumbers.size() == 1){
            set.add(remainingNumbers.get(0));
            superSet.add(set);
        }
        
        for(int i = 0; i < remainingNumbers.size(); i++){
            List<Integer> newSet = new ArrayList<>(set);
            newSet.add(remainingNumbers.get(i));
            
            List<Integer> remainingNums = new ArrayList(remainingNumbers);
            remainingNums.remove(remainingNumbers.get(i));
            
            recursivePermute(remainingNums , superSet, newSet);
        }
    }
    
    public static void reorderList(ListNode head) {
        if(head == null)
            return;
        ListNode pointer = head;
        ListNode doublePointer = head.next;
        while(doublePointer != null){
            pointer = pointer.next;
            doublePointer = doublePointer.next;
            if(doublePointer != null)
                doublePointer = doublePointer.next;
        }
        Stack<ListNode> stack = new Stack<>();
        while(pointer != null){
            stack.push(pointer);
            pointer = pointer.next;
            stack.peek().next = null;
        }
        
        pointer = head;
        while(!stack.isEmpty()){
            ListNode temp = pointer.next;
            pointer.next = stack.pop();
            pointer.next.next = temp;
            pointer = temp;
        }
        pointer.next = null;
    }
    
    public static int[] searchRange(int[] nums, int target) {
        int low = -1;
        int high = nums.length;
        int index = -1;
        while(low + 1 < high){
            int middle = (low + high) / 2;
            if(nums[middle] == target){
                index = middle;
                break;
            }
            if(nums[middle] < target)
                low = middle;
            else
                high = middle;
        }
        
        if(index == -1)
            return new int[]{-1, -1};
        
        //traverseleft
        int leftIndex = index;
        while(leftIndex > 0){
            if(nums[leftIndex - 1] != target){
                break;
            }
            leftIndex++;
        }
        
        int rightIndex = index;
        while(rightIndex < nums.length - 1){
            if(nums[rightIndex + 1] != target)
                break;
            rightIndex++;
        }
        
        return new int[]{leftIndex, rightIndex};
    }
    
    public static double myPow(double x, int n) {
        if(n == 1)
            return x;
        if(n == -1)
            return 1/x;
        if(n == 0)
            return 1;
        int index = 4;
        int absN = Math.abs(n);
        double returnVal = x * x;
        Map<Integer, Double> map = new HashMap<>();
        map.put(2, returnVal);
        while(index <= absN){
            map.put(index, map.get(index / 2) * map.get(index / 2));
            index = index * 2;
        }
        
        index = index / 2;
        returnVal = map.get(index);
        
        while(index != absN){
            returnVal = returnVal * x;
            index++;
        }
        return n < 0 ? 1/returnVal : returnVal;
    }
    
    public static boolean wordBreak(String s, List<String> wordDict) {
        Map<String, Integer> map = new HashMap<>();
        wordDict.forEach(word -> map.put(word, 1));
        return checkWordRecursively(s, map);
    }
    
    private static boolean checkWordRecursively(String s, Map<String, Integer> map) {
        if(s.isEmpty())
            return true;
        Set<String> words = new HashSet<>(map.keySet());
        for (String word : words){
            map.remove(word);
            boolean valid = checkWordRecursively(s.replaceAll(word, "").trim(), map);
            map.put(word, 1);
            if(valid)
                return true; 
        }
        
        return false;
    }
    
    
    public static int numDecodings(String s) {
        int length = s.length();
        Map<Integer, Integer> ways = new HashMap<>();
        return recursion(s, 0 , length, ways);
    }
    
    
    public static int recursion(String s, int index ,int length, Map<Integer, Integer> ways){
        if(index == length){
            return 1;
        }
        
        if(s.charAt(index) == '0')
            return 0;
        
        if(index == length - 1){
            return 1;
        }
        if(ways.containsKey(index))
            return ways.get(index);
        
        int total = 0;
        //current one
        total += recursion(s, index + 1, length, ways);
        String smallString = s.substring(index, index + 2 < length ? index +2 : length);
        int twoChar = 0;
        if(Integer.parseInt(smallString) > 26 || Integer.parseInt(smallString) < 10){
            twoChar = 0;
        }
        else{
            total += recursion(s, index + 2, length, ways);
        }
        
        ways.put(index, total);
        return total;
    }
    
    public static void main(String[] args) {
        Map<String, Node> result = readJson();
        wordBreak("catsandog", Arrays.asList("cats","dog","sand","and","cat") );
        int decodingWays = numDecodings("1112");
        myPow(0.00001 , 2147483647);
        searchRange(new int[]{1}, 1);
        search(new int[]{1,3, 5}, 3);
        divide(100000
                ,10);
        adjacentSum(0, 0, new int[]{2, 4, 6, 2, 5});
        removeInvalidParentheses("()())()");
        hasPath(new long[][]{{0, 1, 1}, {0, 1, 1}, {0, 0, 0}}, 3);
        minimumOccurrences(new long[]{10,941,13,13, 941});
        maxProfit(new long[]{-1, 6, 7 , -2, 10});
        permute(new int[]{1, 2, 3, 4} );
       Map<String, String> map = new HashMap<>();
        ListNode l1 = new ListNode(1);
        ListNode l2 = new ListNode(2);
        ListNode l3 = new ListNode(3);
        ListNode l4 = new ListNode(4);
    
        l1.next = l2;
        l2.next = l3;
        l3.next = l4;
        
        reorderList(l1);
        //mergeTwoLists(l1, l2);
        twoSum(new int[]{2, 3, 4}, 6);
        threeSum(new int[]{-1,0,1,2,-1,-4});
        romanToInt("III");
        myAtoi("4193 with words");
        //longestPalinString("abcdefg");
        //longestPalindrome("aaa");
        //    findMedianSortedArrays(new int[]{1,2}, new int[]{3,4});
        //uniquePathsWithObstacles(new int[][]{{0,0,0},{0,1,0},{0,0,0}});

       // magicLand(new int[]{2,4 , 2,3,1}, 3,5);
        System.out.println(letterCombinations("23"));
    }
}
