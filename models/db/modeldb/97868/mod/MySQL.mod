: $Id: MySQL.mod,v 1.24 2007/04/27 00:41:06 samn Exp $ 
:
: samn at neurosim dot downstate dot edu
:
: if you make fixes, updates, or modifications please send
: me a copy

:thanks to http://www.geocities.com/jahan.geo/mysql_c_by_example.html
:and to the dev.mysql.com docs for MySQL API examples

NEURON {
  SUFFIX mysql
  GLOBAL verbose
}

PARAMETER {
  verbose = 0
}

VERBATIM

#include <stdlib.h>

#include "mysql/mysql.h"
#include "mysql/errmsg.h"

MYSQL g_mysql;
static int g_iInit = 0;
MYSQL_RES* g_result = 0;

char g_user_name[1024]={0};
char g_user_pass[1024]={0};
char g_host[1024]={0};

extern Object** hoc_objgetarg();
extern int ivoc_list_count(Object*);
int list_vector_px();
extern char* hoc_object_name(Object*);

void FreeRes(MYSQL_RES** ppRes){
  if(!ppRes) return;
  mysql_free_result(*ppRes);
  *ppRes = 0;
}

int ConnectReal(){
  mysql_options(&g_mysql,MYSQL_READ_DEFAULT_GROUP,"nrniv");
  //turn on auto-reconnect
  my_bool rt = 1;
  mysql_options(&g_mysql,MYSQL_OPT_RECONNECT , &rt);
  if(!mysql_real_connect(&g_mysql,g_host,g_user_name,g_user_pass,NULL,0,NULL,0)){
    printf("\nConnectReal ERRA: Failed to connect to MySQL server\n");
    return 0.0;
  }
  if(0!=mysql_autocommit(&g_mysql,1)){
    printf("\nConnectReal ERRB: Couldn't turn auto-commit mode on\n");
    return 0.0;
  }
  return 1.0;
}

int InitReal(){
  g_iInit = 0;
  if(mysql_init(&g_mysql)==NULL){
    printf("\nInitReal ERRA: Failed to init MySQL connection\n");
    return 0.0;
  } else {
    g_mysql.reconnect=1;//turn on auto-reconnect
    if(ConnectReal()){
      g_iInit=1;
      return 1.0;
    }
    return 0.0;
  }
}

int CheckConnection(){
  if(!g_iInit){
    printf("MySQL ERRA: never connected to MySQL server\n");
    return 0;
  }  
  g_mysql.reconnect=1;//make sure auto-reconnect on
  int iCheck = mysql_ping(&g_mysql); //check if connection timed-out  
  if(iCheck==0) return 1.0; //server still connected  
  if(iCheck==CR_SERVER_GONE_ERROR){  
    printf("MySQL Connection timed out, reconnecting...\n");
    return 1;
  }
  printf( "Lost connection to MySQL server: Error: %s\n", mysql_error(&g_mysql));
  return 0;
}

ENDVERBATIM

:check reconnection option
FUNCTION ReCon(){
  VERBATIM
  printf("g_mysql.reconnect=%d\n",g_mysql.reconnect);
  return g_mysql.reconnect;
  ENDVERBATIM
}

: closes any open connection to MySQL server
FUNCTION Close(){
  VERBATIM
  if(!g_iInit){
    printf("\nCloseMySQL ERRA: No connection to close\n");
    return 0;
  }
  mysql_close(&g_mysql);
  g_iInit = 0;
  printf("\nClosed MySQL connection\n");
  return 1.0;
  ENDVERBATIM
}


: initialize MySQL engine & connect to MySQL server
: user must supply host-name , user-name, password
: to connect to MySQL server
: returns 1.0 iff successful
: Init(host,user,pass)
FUNCTION Init(){
  VERBATIM
  if(g_iInit){
    printf("Init: Already initialized, checking connection status...\n");
    return (double) CheckConnection();
  }

  if(!ifarg(3)){
    printf("\nInit ERRB: usage: Init(host,user,passwd)\n");
    return 0.0;
  }

  char* tmp  = gargstr(1);
  int iLen = strlen(tmp);
  if(iLen>=1024){
    printf("Init ERR: host string must be < 1024 chars %d!\n",iLen);
    return 0.0;
  } else {
    strcpy(g_host,tmp);
  }

  tmp = gargstr(2);
  iLen = strlen(tmp);
  if(iLen>=1024){
    printf("Init ERR: user string must be < 1024 chars %d!\n",iLen);
    return 0.0;
  } else {
    strcpy(g_user_name,tmp);
  }

  tmp = gargstr(3);
  iLen = strlen(tmp);
  if(iLen>=1024){
    printf("Init ERR: passwd string must be < 1024 chars %d!\n",iLen);
    return 0.0;
  } else {
    strcpy(g_user_pass,tmp);
  }

  return (double) InitReal();

  ENDVERBATIM
}

: SelectDB(dbname)
: returns 1.0 iff successful
FUNCTION SelectDB(){
  VERBATIM
 
  if(!CheckConnection()) return 0.0;

  if(!ifarg(1)){
    printf("SelectDB ERRA: must supply database name\n");
    return 0.0;
  }

  char* dbname = gargstr(1);

  if(mysql_select_db(&g_mysql,dbname)==0){
    printf( "Database %s Selected\n",dbname);
  } else {
    printf( "Failed to connect to Database: Error: %s\n", mysql_error(&g_mysql));
  }

  return 1.0;
  
  ENDVERBATIM
}

: frees results of Select, responsibility of hoc user
FUNCTION FreeResults(){
  VERBATIM

  if(!CheckConnection()) return 0.0;

  if(!g_result){
    printf("FreeResults ERRA: No results to free\n");
    return 0.0;
  }

  FreeRes(&g_result);

  return 1.0;

  ENDVERBATIM
}

: check # of columns from previous Select call
FUNCTION NumCols(){
  VERBATIM

  if(!CheckConnection()) return 0.0;

  if(g_result){
    return (double) mysql_num_fields(g_result);
  } else {
    printf("NumCols ERRA: no result set\n");
    return 0.0;
  }

  ENDVERBATIM
}

: check # of rows from previous Select call
FUNCTION NumRows(){
  VERBATIM

  if(!CheckConnection()) return 0.0;

  if(g_result){
    return (double) mysql_num_rows(g_result);
  } else {
    printf("NumRows ERRA: no result set\n");
    return 0.0;
  }

  ENDVERBATIM
}

: get rows from previous Select
: into list of vectors (each vec is dimension/column)
: returns -1.0 on error, otherwise number of rows
FUNCTION GetRows(){
  VERBATIM

  if(!CheckConnection()) return 0.0;

  if(!ifarg(1)){
    printf("GetRows ERRA: must pass in list of vectors as arg 1\n");
    return -1.0;
  }

  if(!g_result){
    printf("GetRows ERRB: no SQL results\n");
    return -1.0;
  }

  unsigned int num_fields = mysql_num_fields(g_result);
  unsigned int num_rows = mysql_num_rows(g_result);

  if(num_fields == 0){
    printf("GetRows ERRC: empty SQL results rows=%d fields=%d\n",num_rows,num_fields);
    return -1.0;
  }

  Object* pList = *hoc_objgetarg(1);
  if(!pList){
    printf("GetRows ERRD: first arg must be list of vectors\n");
    return -1.0;
  }

  int iColumns = ivoc_list_count(pList);

  if(iColumns < num_fields){
    printf("GetRows ERRE: num columns in list < num SQL select columns\n");
    return -1.0;
  }

  int i;
  double** vvo = (double**) malloc(iColumns * sizeof(double*));
  if(!vvo){
    printf("GetRows ERRF: out of memory\n");
    return -1.0;
  }

  int iVecSz = 0;
  for(i=0;i<iColumns;i++){
    iVecSz = list_vector_px(pList,i,&vvo[i]);
    if(iVecSz < num_rows){
      printf("GetRows ERRG: vec size %d < SQL result rows %d\n",iVecSz,num_rows);
      free(vvo);
      return -1.0;
    }
  }

  // retrieve rows
  MYSQL_ROW row; int j = 0;
  while ((row = mysql_fetch_row(g_result))) { 
    unsigned long* lengths = mysql_fetch_lengths(g_result);
    int i;
    for(i = 0; i < num_fields; i++) { 
      vvo[i][j] = atof(row[i]);
      if(verbose) printf("[%.*s] ", (int) lengths[i], row[i] ? row[i] : "NULL"); 
    } 
    if(verbose) printf("\n"); 
    j++;
  }

  free(vvo);

  return (double) num_rows;  

  ENDVERBATIM
}

VERBATIM

int StringEqualIgnoreCase(char* p1,char* p2,int n){
  int i;
  for(i=0;i<n;i++){
    if(tolower(p1[i])!=tolower(p2[i])){
      return 0;
    }
  }
  return 1;
}

ENDVERBATIM

: takes vector and returns # of times it exists as a row in table_name
: Find(table_name,Vector) also allows partial row match on first 
: min(vector.size,table.columns) columns stores results in g_result for later
: retrieval returns -1.0 on error, otherwise num_rows found matching vector
FUNCTION Find(){
  VERBATIM

  if(!CheckConnection()) return 0.0;

  if(!ifarg(1)){
    printf("Find ERRA: must supply table name\n");
    return -1.0;
  }

  char* table = gargstr(1);

  if(g_result){
    FreeRes(&g_result);
    if(verbose) printf("Find Warning: previous SQL results freed\n");
  }

  double* pVec;
  int iVecSz = vector_arg_px(2,&pVec);
  if(!pVec || iVecSz < 1){
    printf("Find ERRB: arg 2 must be a Vector > 0 size\n");
    return -1.0;
  }

  MYSQL_RES* res = mysql_list_fields(&g_mysql,table,NULL);
  if(!res){
    printf( "Failed to search for record: Error: %s\n", mysql_error(&g_mysql));
    return -1.0;
  }

  unsigned int num_fields = mysql_num_fields(res);

  if(iVecSz > num_fields){
    printf("Find ERRC: incorrect sized vector %d , %s num cols = %d\n",iVecSz,table,num_fields);
    FreeRes(&res);
    return -1.0;
  }

  char sql_select[2048]={0};
  sprintf(sql_select,"select * from %s where ",table);
  char sql_val[1024]={0};

  int i; 

  int iSearchCols = num_fields < iVecSz ? num_fields : iVecSz;

  for(i = 0; i < num_fields && i < iVecSz; i++) { 
    MYSQL_FIELD* field = mysql_fetch_field(res);
    sprintf(sql_val,"%s=%f",field->name,pVec[i]);
    strcat(sql_select,sql_val);
    if(i < iSearchCols - 1 ){
      strcat(sql_select," and ");
    } else {
      strcat(sql_select,";");
    }
  } 

  FreeRes(&res);

  if(verbose) printf("Find: sql command = %s\n",sql_select);

  if(0==mysql_real_query(&g_mysql,sql_select,strlen(sql_select))){
    g_result = mysql_store_result(&g_mysql); //store whole row for later retrieval
    if(g_result){
      unsigned int num_rows = mysql_num_rows(g_result);
      return (double) num_rows;
    } else {
      fprintf(stderr,"Find ERRE: %s\n",mysql_error(&g_mysql));
      return -1.0;
    }
  } else {
    printf("Find ERRD: %s\n",mysql_error(&g_mysql));
    return -1.0;
  }
 

  ENDVERBATIM
}


: updates a single col of
: a table. 
: UpdateCol(table_name,col_name,order_by_column_name,vector_of_values,start_idx)
: col_name is the column that will be updated
: order_by_column_name is the column that stores
: ids, if no column stores ids, updating a column
: does not make much sense because the storage
: order of column values may not be what the
: user is expecting. start_idx is starting value
: of order by column index. it is incremented
: for each row of a column.
: so it will be 
: update table set col_name = vec[0] where order_by_column_name=start_idx;
: update table set col_name = vec[1] where order_by_column_name=start_idx+1;
:   ...
: update table set col_name = vec[n] where order_by_column_name=start_idx+n;
FUNCTION UpdateCol(){
  VERBATIM

  if(!CheckConnection()) return 0.0;

  if(!ifarg(4)){
    printf("UpdateCol ERRA: usage UpdateCol(table_name,col_name,order_by_col_name,vec_of_vals\n");
    return 0.0;
  }

  char* table = gargstr(1);
  char* col_name = gargstr(2);
  char* order_by_col = gargstr(3);

  double* pVec = 0;
  int iSz = vector_arg_px(4,&pVec);
  if(!pVec){
    printf("UpdateCol ERRB: couldn't get vector arg 4!\n");
    return 0.0;
  }

  int idx = ifarg(5) ? (int) *getarg(5) : 0;

  char sql_command[8192]={0};

  int i;
  for(i=0;i<iSz;i++){

    sprintf(sql_command,"update %s set %s=%f where %s=%d;",table,col_name,pVec[i],order_by_col,idx);

    idx++;

    if(verbose) printf("sql cmd = %s\n",sql_command);

    if(mysql_real_query(&g_mysql,sql_command,strlen(sql_command))==0){
      if(verbose) printf( "col %s updated\n");
    } else {
      printf( "Failed to update record: Error: %s\n", mysql_error(&g_mysql));
    }
  }

  return 1.0;
  ENDVERBATIM
}

: inserts data into existing table
: Insert(table_name,list_of_vectors or vector)
: Vector should have same size as # of columns in table , so Insert
: will add 1 row for Vector arg if arg is List, it should have
: num_cols vectors and vectors.size rows will be inserted into table
: returns 1.0 iff success
FUNCTION Insert(){

  VERBATIM

  if(!CheckConnection()) return 0.0;

  if(!ifarg(1)){
    printf("Insert ERRB: must supply table name\n");
    return 0.0;
  }

  char* table = gargstr(1);

  if(g_result){
    FreeRes(&g_result);
    if(verbose) printf("Insert Warning: previous SQL results freed\n");
  }

  Object* pObj = *hoc_objgetarg(2);
  if(!pObj){
    printf("Insert ERRC: second arg must be List of Vector(s) or Vector\n");
    return 0.0;
  }
 
  char sql_insert[2048] = {0};
  sprintf(sql_insert,"insert into %s values( ",table);
 
  char sql_command[8192] = {0};  
  char sql_val[256] = {0};
  
  if(!strncmp(hoc_object_name(pObj),"List",4)){
 
    int iColumns = ivoc_list_count(pObj);
 
    if(iColumns == 0){
      printf("Insert ERRD: empty list\n!");
      return 0.0;
    }
 
    int i;
    double** vvo = (double**) malloc(iColumns * sizeof(double*));
    if(!vvo){
      printf("Insert ERRE: out of memory\n");
      return 0.0;
    }
 
    int iVecSz = list_vector_px(pObj,0,&vvo[0]);
    for(i=1;i<iColumns;i++){
      int iTmp = list_vector_px(pObj,i,&vvo[i]);
      if(iTmp != iVecSz){
        printf("Insert ERRF: vectors of different sizes %d %d!\n",iVecSz,iTmp);
        free(vvo);
        return 0.0;
      }
    }
 
    int j = 0;
 
    for(i=0;i<iVecSz;i++){
      strcpy(sql_command,sql_insert);
      for(j=0;j<iColumns;j++){
        sprintf(sql_val,"%f",vvo[j][i]);
        strcat(sql_command,sql_val);
        if(j<iColumns-1) strcat(sql_command,",");
      }
      strcat(sql_command," );");
 
      if(verbose) printf("sql_command = %s\n",sql_command);
 
      if(mysql_real_query(&g_mysql,sql_command,strlen(sql_command))==0){
        if(verbose) printf( "Record Added\n");
      } else {
        printf( "Failed to add record: Error: %s\n", mysql_error(&g_mysql));
      }
    }
 
    free(vvo);
  } else if(!strncmp(hoc_object_name(pObj),"Vector",6)){
    double* pVec = NULL;
    int iCols = vector_arg_px(2,&pVec);
    int i = 0;
    strcpy(sql_command,sql_insert);
    for(i=0;i<iCols;i++){
      sprintf(sql_val,"%f",pVec[i]);
      strcat(sql_command,sql_val);
      if(i<iCols-1) strcat(sql_command,",");
    }
    strcat(sql_command," );");
    if(verbose) printf("sql_command = %s\n",sql_command);

    if(mysql_real_query(&g_mysql,sql_command,strlen(sql_command))==0){
      if(verbose) printf("Record Added\n");
    } else {
      printf("Failed to add record: Error: %s\n",mysql_error(&g_mysql));
    }
  } else {
    printf("Insert ERRG: Invalid object type %s\n",hoc_object_name(pObj));
    return 0.0;
  }

  return 1.0;
  
  ENDVERBATIM
}

: does a sql select and keeps results around for hoc user to retrieve
: hoc user must free results at a certain point
: returns -1.0 on error, otherwise # of rows found
FUNCTION Select(){

  VERBATIM

  if(!CheckConnection()) return 0.0;

  if(!ifarg(1)){
    printf("SelectSQL ERRA: must supply SQL select string!\n");
    return -1.0;
  }

  char* query = gargstr(1);

  if(strlen(query)<6 || !StringEqualIgnoreCase(query,"select",6)){
    printf("SelectSQL ERRB: must supply SQL select string!\n");
    return -1.0;
  }

  //must free previous results
  if(g_result){
    FreeRes(&g_result);
    if(verbose) printf("SelectSQL Warning: previous SQL results freed\n");
  }  

  if(0==mysql_real_query(&g_mysql,query,strlen(query))){

    g_result = mysql_store_result(&g_mysql);
 
    // there are rows
    if (g_result)  {

      unsigned int num_fields = mysql_num_fields(g_result);
      unsigned int num_rows = mysql_num_rows(g_result);

      printf("%d rows, %d fields per row, call GetRows to get\n",num_rows,num_fields);

      return (double) num_rows;

    } else {
        // mysql_store_result() should have returned data
        fprintf(stderr, "SelectSQL ERRC : %s\n", mysql_error(&g_mysql));
        return -1.0;
      }
  } else {
    printf("SelectSQL ERRB: %s\n",mysql_error(&g_mysql));
    return -1.0;
  }

  ENDVERBATIM
}

: gets col names from last sql select must pass in correct # of
: char* 's and they must have sufficient length to store col names
FUNCTION GetColNames(){
  VERBATIM

  if(!CheckConnection()) return 0.0;

  if(!g_result){
    printf("GetColNames ERRA: no sql select results\n"); 
    return 0.0;  
  }

  unsigned int num_fields = mysql_num_fields(g_result);

  if(!ifarg(num_fields)){ 
    printf("GetColNames ERRB: must pass in at least %d char*'s\n",num_fields);
    return 0.0;
  }

  MYSQL_FIELD* fields = mysql_fetch_fields(g_result);
  if(!fields){
    printf("GetColNames ERRC: couldn't get column names\n");
    return 0.0;
  }

  int i;
  for(i=0;i<num_fields;i++){
    strcpy(gargstr(i+1),fields[i].name);
  }

  return 1.0;

  ENDVERBATIM
}

VERBATIM

double PrintRows(MYSQL_RES* res){

  if(!res) return 0.0;

  unsigned int num_fields = mysql_num_fields(res);
  unsigned int num_rows = mysql_num_rows(res);

  MYSQL_FIELD* fields = mysql_fetch_fields(res);

  // print column names
  printf("\n__________________________________________________\n");
  int i;
  for(i=0;i<num_fields;i++){
    printf("[%s]\t",fields[i].name);
  }
  printf("\n__________________________________________________\n\n");

  // print rows
  MYSQL_ROW row;
  while ((row = mysql_fetch_row(res))) { 
    unsigned long* lengths = mysql_fetch_lengths(res);
    int i;
    for(i = 0; i < num_fields; i++) { 
      printf("[%.*s]\t", (int) lengths[i], row[i] ? row[i] : "NULL"); 
    } 
    printf("\n"); 
  }     
  printf("\n__________________________________________________\n");
  printf("%d row(s), %d col(s)\n",num_rows,num_fields);
  printf("__________________________________________________\n\n");

  return (double) num_rows;
}

ENDVERBATIM

: lists all available dbs
: returns -1 iff error, otherwise # of dbs
FUNCTION ListDBs(){
  VERBATIM

  if(!CheckConnection()) return 0.0;

  if(g_result){
    FreeRes(&g_result);
    if(verbose) printf("ListDBs Warning: previous SQL results freed\n");
  }

  MYSQL_RES* dbres = mysql_list_dbs(&g_mysql,NULL);

  if(dbres){

    double num_rows = PrintRows(dbres);    

    FreeRes(&dbres);

    return num_rows;

  } else {

    fprintf(stderr, "ListDBs ERRA : %s\n", mysql_error(&g_mysql));
    return -1.0;

  }
  

  ENDVERBATIM
}

: executes a sql command but doesnt store results for hoc user
: can execute any type of SQL command, i.e. create,select,insert,etc.
: returns -1.0 on error
: Query(query_string)
FUNCTION Query(){
  VERBATIM

  if(!CheckConnection()) return 0.0;

  if(!ifarg(1)){
    printf("Query ERRA: must supply SQL query string!\n");
    return -1.0;
  }

  char* query = gargstr(1);

  //must free previous results
  if(g_result){
    FreeRes(&g_result);
    if(verbose) printf("Query Warning: previous SQL results freed\n");
  }

  if(0==mysql_real_query(&g_mysql,query,strlen(query)))
  {
    MYSQL_RES* result = mysql_store_result(&g_mysql);

    // there are rows
    if(result)
    {    
      double num_rows = PrintRows(result);

      FreeRes(&result);

      return num_rows;
    }
    else
    {
      //mysql_store_result() returned nothing; should it have?
      if(mysql_field_count(&g_mysql) == 0)
      {
        // query does not return data
        // (it was not a SELECT)
        unsigned int num_rows = mysql_affected_rows(&g_mysql);
        printf("Query: affected rows = %d\n",num_rows);
        return (double) num_rows;
      }
      else
      {
        //mysql_store_result() should have returned data
        fprintf(stderr, "Query ERRC : %s\n", mysql_error(&g_mysql));
        return -1.0;
      }
    }    
  }
  else
  {
    printf("Query ERRB: %s\n",mysql_error(&g_mysql));
    return -1.0;
  }

  ENDVERBATIM
}

: display client & server versions
PROCEDURE VersionInfo(){
  VERBATIM
  if(!CheckConnection()) return;
  printf("MySQL Client Version is %s\n",mysql_get_client_info());
  printf("MySQL Server Version is %s\n",mysql_get_server_info(&g_mysql));
  ENDVERBATIM
}


