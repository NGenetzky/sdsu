/*
 * File:   main.c
 * Author: Nathansen
 *
 * PORTA is used for Anlog Sensor
 *  A3 Thermister Reading
 * PORTD is used for LED.
 *  D0  Green LED
 *  D1  Amber LED
 *  D2  Red LED
 *  D3  Run LED
 *  D4  Heater LED
 * PORTE is used for keypad.
 *  E0  COL0
 *  E1  COL1
 *  E2  COL2
 *  E3  COL3
 *  E4  ROW0
 *  E5  ROW1
 *  E6  ROW2
 *  E7  ROW3
 * PORTH
 *  J0  Run Switch
 *  J1  Input Switch
 * Created on January 24, 2015, 3:13 PM
 */
//BOIL
//SIMMER
//SOUP      T=90C, PIC=85C, ADC=F70
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <p18f8722.h>
//C includes
#include <limits.h>
#include <stdint.h>
#include <string.h>
//Device Includes
#include <xc.h>  //Device access
#include <spi.h>
#include <delays.h>
//define for PIC18
#define _XTAL_FREQ 10000000
// configuration bits
#pragma config OSC = HSPLL    // Oscillator Selection bits (HS oscillator)
#pragma config FCMEN = OFF // Fail-Safe Clock Monitor Enable bit (Fail-Safe Clock Monitor disabled)
#pragma config WDT = OFF   // Watchdog Timer (WDT disabled (control is placed on the SWDTEN bit))
#pragma config LVP = OFF   //Single-Supply ISCP Enable bit

//Include personal headers
#include "LCD_Header.h"

//Main Functions
void interrupt high_priority interrupt_at_high_vector(void);
void interrupt low_priority interrupt_at_low_vector(void);
void Setup();
void loop_feedback();
void loop_setValue();
//<<    LCD Functions
#define LCD_POS_LINE2 40
void LCD_updateDisp(int disp);
void LCD_printUInt(int);
void LCD_printDouble(double);
void LCD_disp_Buffer();
void LCD_disp_A3andVoltage();
void LCD_disp_SVandPV();
int LCD_pos;
#define LCD_DISP_INPUT 0
#define LCD_DISP_OUTPUT 1
//      LCD Functions       >>

//<< Tea Feedback Control Project
#define TEA_SETTRIS TRISJ = 0b000000011;
#define PORT_LED_GREEN PORTJbits.RJ3
#define PORT_LED_AMBER PORTJbits.RJ4
#define PORT_LED_RED PORTJbits.RJ5
#define PORT_LED_RUN PORTJbits.RJ6
#define PORT_HEATER PORTJbits.RJ7
#define TEA_ARRAYSIZE 20
#define PORT_SWITCH_RUN PORTJbits.RJ0
#define PORT_SWITCH_INPUT PORTJbits.RJ1
#define STATE_INPUT 0
#define STATE_OUTPUT 1
void TEA_checkTea();
void TEA_detectState();
void TEA_checkTemp();
double TEA_HighTemp, TEA_LowTemp;
double TEA_temperature;
double TEA_RecentTemps[TEA_ARRAYSIZE];
int TEA_TempArrayIndex;
int TEA_State;
//   Tea Feedback Control Project   >>

//<<    Buffer Functions
#define BUFFER_SIZE 20
void    BUFFER_add(char key);
void    BUFFER_clear();
int     BUFFER_toInt();
char    *BUFFER_ptr = NULL;
char    BUFFER_array[BUFFER_SIZE];
//      Buffer Funtions      >>

//<<    Keypad DEFINES and Function declares
#define PORT_KEYPAD PORTE
#define TRIS_KEYPAD TRISE
#define PORT_KEYPAD_C0 PORTEbits.RE0
#define PORT_KEYPAD_C1 PORTEbits.RE1
#define PORT_KEYPAD_C2 PORTEbits.RE2
#define PORT_KEYPAD_C3 PORTEbits.RE3
#define PORT_KEYPAD_R0 PORTEbits.RE4
#define PORT_KEYPAD_R1 PORTEbits.RE5
#define PORT_KEYPAD_R2 PORTEbits.RE6
#define PORT_KEYPAD_R3 PORTEbits.RE7
#define NOKEY 255
void    KEYPAD_Setup();
void    KEYPAD_Loop();
void    KEYPAD_processKeyPress();
char    KEYPAD_getKeyBytewise();
char    KEYPAD_getKey();
char    KEYPAD_getColPressed();
void    KEYPAD_delay();
char    KEYPAD_IdentifyKey(char row, char col);
char    KEYPAD_KeyToASCI(char id);
char    KEYPAD_key2ASCI[] ={'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F'};
char    KEYPAD_lastPressed;
char    KEYPAD_drivenValue;
//  Keypad DEFINES and Function declares     >>

//<<    A2D Function Declares
void    A2D_Loop();
void    A2D_ReadA3ToPortD();
unsigned int     A2D_readResult();
void    A2D_Setup();
void    A2D_start();
int     A2D_getBothByteResult();
int     A2D_getHighByteResult();
double  A2D_getVoltFromDigital(char highByte);
double  A2D_getTemperatureFromVolt(double voltage);
double  A2D_getTempFromA2D();
double  A2D_getVoltFrom8bitA2D();
double  A2D_getVoltFrom10bitA2D();
char    A2D_ResultText[10];
//      A2D Function Declares       >>
//<<        Delay functions
void delayBetweenMeasure();
void delay500ms();
void delayXSec(int loopNum);
//        Delay functions     >>
//ISR
void interrupt high_priority interrupt_at_high_vector(void){

}
void interrupt low_priority interrupt_at_low_vector(void){

}
void Setup(){
    TRISD = 0x00; //Output
    PORTD = 0x00; //Clear to 0
    TEA_SETTRIS //TRISJ = 0xFF; //Input
    TRISH = 0xFF; //Input
    KEYPAD_Setup();
    lcd();
    LCD_pos=0;
    BUFFER_clear();
    A2D_Setup();
    TEA_TempArrayIndex=0;

    BUFFER_add('2');BUFFER_add('0');BUFFER_add('\0');
    PORT_HEATER=0;
    PORT_LED_RED =0;
    PORT_LED_AMBER =0;
    PORT_LED_GREEN =0;
}
int main(int argc, char** argv) {
    Setup();
    //PORT_LED_RUN =1; //Running
    while(1){//LOOP
        //loop_feedback();
        //loop_setValue();
        TEA_detectState();
        if(!PORT_SWITCH_RUN){
            PORT_LED_RUN = 0;
        }
        else if(TEA_State==STATE_INPUT){
            PORT_LED_RUN = 1;
            PORT_HEATER=0;
            PORT_LED_RED =0;
            PORT_LED_AMBER =0;
            PORT_LED_GREEN =0;
            loop_setValue();
        }else if(TEA_State==STATE_OUTPUT){
            PORT_LED_RUN = 1;
            loop_feedback();
        }
    }
}
void loop_setValue(){
    LCD_updateDisp(LCD_DISP_INPUT);
    KEYPAD_processKeyPress();
}
void loop_feedback(){
    TEA_LowTemp = BUFFER_toInt();
    TEA_HighTemp = TEA_LowTemp+2;
    TEA_checkTea();
        //A2D_getHighByteResult();
    LCD_updateDisp(LCD_DISP_OUTPUT);
    delayBetweenMeasure();//delay500ms();
}
void LCD_updateDisp(int disp){
    lcd();
    switch(disp){
        case LCD_DISP_INPUT:
            lcdGoTo(0);
                LCD_disp_Buffer();
            lcdGoTo(40); //Line2
                LCD_disp_SVandPV();
            lcdGoTo(LCD_pos);
            break;
        case LCD_DISP_OUTPUT:
            lcdGoTo(0);
                LCD_disp_A3andVoltage();//LCD_disp_Buffer();
            lcdGoTo(40); //Line2
                LCD_disp_SVandPV();
            lcdGoTo(0);
            break;
    }
}
void LCD_disp_Buffer(){
        for(int i=0; (i<17)&&(i<BUFFER_SIZE);i++){
            lcdChar(BUFFER_array[i]);
        }
}
void LCD_disp_A3andVoltage(){
    lcdWriteString("A3(");
    //LCD_printDouble(A2D_getVoltFrom8bitA2D());
    LCD_printUInt(ADRESH);
    lcdWriteString(") ");
    LCD_printDouble(A2D_getVoltFrom8bitA2D());
    lcdChar('V');
}
void LCD_disp_SVandPV(){
    lcdWriteString("SV:");
    LCD_printUInt(BUFFER_toInt());
    //lcdChar(')');
    lcdWriteString(" PV:");
    LCD_printDouble(TEA_temperature);
    //lcdChar('');
}
/*
    //lcdWriteString("V(");
    lcdGoTo(40);
    lcdWriteString("Key(");
    lcdChar(KEYPAD_KeyToASCI(KEYPAD_lastPressed));
    lcdChar(')');
    lcdGoTo(47); //Line2
    lcdWriteString("AN3(");
    lcdWriteString(A2D_ResultText);
    lcdChar(')');
    */
void LCD_printUInt(int intNum){
    int counter;
    unsigned char num[8];
    sprintf(num, "%u", intNum);
    for(counter = 0; counter < strlen(num); counter++)
    {
        lcdChar(num[counter]);
    }
}
void LCD_printDouble(double doubleNum){
    int counter;
    unsigned char num[10];
    sprintf(num, "%.2f", doubleNum);
    for(counter = 0; counter < strlen(num); counter++)
    {
        lcdChar(num[counter]);
    }
}
//<<   Tea Feedback Control Project
void TEA_checkTea(){
    TEA_checkTemp();
    if(TEA_temperature<TEA_LowTemp-2){
        PORT_HEATER=1;
        PORT_LED_RED =0;
        PORT_LED_AMBER =1;
        PORT_LED_GREEN =0;
    }
    else if((TEA_LowTemp-2<TEA_temperature)&&(TEA_temperature<TEA_LowTemp)){
        //Partially on
        if(TEA_temperature<90){
        PORT_HEATER=1;
        __delay_ms(50);
        }

        PORT_LED_RED =0;
        PORT_LED_AMBER =1;
        PORT_LED_GREEN =0;
        PORT_HEATER=0;
    }
    else if((TEA_LowTemp<TEA_temperature)&&(TEA_temperature<TEA_HighTemp)){
        //Partially on
        PORT_HEATER=1;
        __delay_ms(25);

        PORT_HEATER=0;
        PORT_LED_RED =0;
        PORT_LED_AMBER =0;
        PORT_LED_GREEN =1;
    }
    else if(TEA_HighTemp<TEA_temperature){
        PORT_HEATER=0;
        PORT_LED_RED =1;
        PORT_LED_AMBER =0;
        PORT_LED_GREEN =0;
    }
}
void TEA_checkTemp(){
    A2D_getHighByteResult();
    TEA_RecentTemps[TEA_TempArrayIndex] = A2D_getTempFromA2D();
    if(TEA_TempArrayIndex<TEA_ARRAYSIZE){
        TEA_TempArrayIndex++;
    }else{
        TEA_TempArrayIndex =0;
    }
    double sum= 0;
    for(int i=0; i< TEA_ARRAYSIZE ; i++){
        sum = sum + TEA_RecentTemps[i];
    }
    TEA_temperature = sum/TEA_ARRAYSIZE;
}
void TEA_detectState(){
    if(PORT_SWITCH_INPUT){
            KEYPAD_delay();
            if(PORT_SWITCH_INPUT){//Debounced
                TEA_State = STATE_INPUT;
            }
        }
        else{
            KEYPAD_delay();
            if(!PORT_SWITCH_INPUT){//Debounced
                TEA_State = STATE_OUTPUT;
            }
        }
}
//   Tea Feedback Control Project   >>
//<<    KEYPAD functions
void KEYPAD_Setup(){
    PORT_KEYPAD = 0x0F; // Force Default, COL=1, ROW =0
    TRIS_KEYPAD = 0x0F; // Low nibble(COL) is input, High nibble(ROW) is output
    TRISD = 0x00;
    //PORTD = 0xAA;
    BUFFER_ptr = BUFFER_array;
}
void KEYPAD_Loop(){
    char mKeypad_key=NOKEY;
    while(1){
        mKeypad_key = KEYPAD_getKey();
        //PORTD = mKeypad_key;
        mKeypad_key = NOKEY;
    }
}

void KEYPAD_processKeyPress(){
    KEYPAD_getKey();
    if(KEYPAD_lastPressed==NOKEY){
        return;
    }
    char charToAdd;
    char numEntryTable[] =
           {'1','2','3',' ',
            '4','5','6',' ',
            '7','8','9',' ',
            ' ','0','.','\0'};
    switch(KEYPAD_lastPressed){
        case 0:
        case 1:
        case 2:
        case 4:
        case 5:
        case 6:
        case 8:
        case 9:
        case 10:
        case 13:
        //case 14: //Decimal Point
            // Regular characters
            charToAdd = numEntryTable[KEYPAD_lastPressed];
            BUFFER_add(charToAdd);
            break;
        case 3://UP
            //lcdCommand(0b00011111); //Shift disp right
            //BUFFER_ptr--;
            //LCD_pos--;
            break;
        case 7://DOWN
            //lcdCommand(0b00011011);//Shift disp left
            //BUFFER_ptr++;
            //LCD_pos++;
            break;
        case 11://2nd
            lcdCommand(0b00001111);
            break;
        case 12://Clear
            BUFFER_clear();
            break;
        case 15:
            charToAdd = '\0';
            BUFFER_add(charToAdd);
            break;
        default:
            charToAdd = NOKEY;// Special functions
    }
    delay500ms();
}
char KEYPAD_poll(){
    PORT_KEYPAD = 0x0F; //Ground all rows
    TRIS_KEYPAD = 0x0F; // Low nibble(COL) is input, High nibble(ROW) is output
    return KEYPAD_getColPressed()!=NOKEY;
}
char KEYPAD_readKey(){
    char colPressed, keyPressed, rowPressed;
    colPressed = NOKEY; keyPressed =NOKEY; rowPressed=NOKEY;
    PORT_KEYPAD = 0xFF; // Pull all rows high.
//Check row by grounding ROW0
    PORT_KEYPAD_R0 = 0; KEYPAD_drivenValue = 0b11101111;
    KEYPAD_delay();
    colPressed = KEYPAD_getColPressed();
    if(colPressed!=255){
        rowPressed = 0;
        KEYPAD_lastPressed = KEYPAD_IdentifyKey(rowPressed, colPressed);
        return KEYPAD_lastPressed;
    }
    PORT_KEYPAD_R0 = 1;
//Check row by grounding ROW1
    PORT_KEYPAD_R1 = 0; KEYPAD_drivenValue = 0b11011111;
    KEYPAD_delay();
    colPressed = KEYPAD_getColPressed();
    if(colPressed!=255){
        rowPressed = 1;
        KEYPAD_lastPressed = KEYPAD_IdentifyKey(rowPressed, colPressed);
        return KEYPAD_lastPressed;
    }
    PORT_KEYPAD_R1 = 1;
//Check row by grounding ROW2
    PORT_KEYPAD_R2 = 0; KEYPAD_drivenValue = 0b10111111;
    KEYPAD_delay();
    colPressed = KEYPAD_getColPressed();
    if(colPressed!=255){
        rowPressed = 2;
        KEYPAD_lastPressed = KEYPAD_IdentifyKey(rowPressed, colPressed);
        return KEYPAD_lastPressed;
    }
    PORT_KEYPAD_R2 = 1;
//Check row by grounding ROW3
    PORT_KEYPAD_R3 = 0; KEYPAD_drivenValue = 0b01111111;
    KEYPAD_delay();
    colPressed = KEYPAD_getColPressed();
    if(colPressed!=255){
        rowPressed = 3;
        KEYPAD_lastPressed = KEYPAD_IdentifyKey(rowPressed, colPressed);
        return KEYPAD_lastPressed;
    }
    PORT_KEYPAD_R3 = 1;
    return NOKEY;
}
char KEYPAD_getKey(){
//    if(KEYPAD_poll()){
//        KEYPAD_lastPressed =KEYPAD_readKey();
//        return KEYPAD_lastPressed;
//    }
//    delay500ms();
    PORT_KEYPAD = 0x0F; //Ground all rows
    TRIS_KEYPAD = 0x0F; // Low nibble(COL) is input, High nibble(ROW) is output
    char colPressed, keyPressed, rowPressed;
    colPressed = NOKEY; keyPressed =NOKEY; rowPressed=NOKEY;
    while((keyPressed==NOKEY)&&(TEA_State==STATE_INPUT)){
        PORT_KEYPAD = 0x0F; //Ground all rows
        TEA_detectState();
        colPressed = KEYPAD_getColPressed();
        if(colPressed==NOKEY){//No Key pressed
            //Check again
        }
        else{ //Key Pressed
            colPressed =NOKEY;
            PORT_KEYPAD = 0xFF; // Pull all rows high.
        //Check row by grounding ROW0
            PORT_KEYPAD_R0 = 0; KEYPAD_drivenValue = 0b11101111;
            KEYPAD_delay();
            colPressed = KEYPAD_getColPressed();
            if(colPressed!=255){
                rowPressed = 0;
                break;
            }
            PORT_KEYPAD_R0 = 1;
        //Check row by grounding ROW1
            PORT_KEYPAD_R1 = 0; KEYPAD_drivenValue = 0b11011111;
            KEYPAD_delay();
            colPressed = KEYPAD_getColPressed();
            if(colPressed!=255){
                rowPressed = 1;
                break;
            }
            PORT_KEYPAD_R1 = 1;
        //Check row by grounding ROW2
            PORT_KEYPAD_R2 = 0; KEYPAD_drivenValue = 0b10111111;
            KEYPAD_delay();
            colPressed = KEYPAD_getColPressed();
            if(colPressed!=255){
                rowPressed = 2;
                break;
            }
            PORT_KEYPAD_R2 = 1;
        //Check row by grounding ROW3
            PORT_KEYPAD_R3 = 0; KEYPAD_drivenValue = 0b01111111;
            KEYPAD_delay();
            colPressed = KEYPAD_getColPressed();
            if(colPressed!=255){
                rowPressed = 3;
                break;
            }
            PORT_KEYPAD_R3 = 1;
        }//End Key Pressed
    }// End While loop
    //The Row has been determined
    KEYPAD_lastPressed = KEYPAD_IdentifyKey(rowPressed, colPressed);
    return KEYPAD_lastPressed;
}
char KEYPAD_getColPressed(){
    char colPressed =255;
    if(!PORT_KEYPAD_C0){
        colPressed = 0;
        }
    else if(!PORT_KEYPAD_C1){
        colPressed = 1;
    }
    else if(!PORT_KEYPAD_C2){
        colPressed = 2;
    }
    else if(!PORT_KEYPAD_C3){
        colPressed = 3;
    }
    else{
        colPressed = 255;
    }
    return colPressed;
}
void KEYPAD_delay(){
    unsigned char i;
    for(i = 0x00; i< 0xFF; i++)
    return;
}
char KEYPAD_IdentifyKey(char rowPressed, char colPressed){
    // 0 1 2 3
    // 4 5 6 7
    // 8 9 A B
    // C D E F
    return rowPressed*4+colPressed;//Give hex value
    //return KEYPAD_KeyToASCI(rowPressed*4+colPressed);//Give ASCI value
}
char KEYPAD_KeyToASCI(char id){
    if(id == NOKEY){return NOKEY;}
    else{return KEYPAD_key2ASCI[id];}
}
//    KEYPAD functions      >>

void BUFFER_add(char key){
    *BUFFER_ptr = key;
    if(BUFFER_ptr>=&BUFFER_array[BUFFER_SIZE-1]){
        BUFFER_ptr = BUFFER_array; //Go back to begining
        LCD_pos =0;
    }
    else{
        BUFFER_ptr++;
        LCD_pos++;
    }
}
void    BUFFER_clear(){
    for(int i=0;i<BUFFER_SIZE;i++){
        BUFFER_array[i]=' ';
    }
    BUFFER_ptr = BUFFER_array;
    LCD_pos=0;
}
int BUFFER_toInt(){
    int iNum=0;
    unsigned char intNum[10];
    char c;
    for(int i=0;i<BUFFER_SIZE;i++){
        c =BUFFER_array[i];
        if( ('0'<=c)&&(c<='9') ){
            intNum[iNum] = c;
            iNum++;
        }
        else if(c=='\0'){
            break; // Don't care about rest of buffer
        }
        else if(c=='.'){
            //Nothing
        }
        else{
            //PORTD = 0x50;
        }
    }
    return strtol (intNum,NULL,10);
}
//<<        A2D functions
void A2D_Loop(){
    TRISD = 0x00;
    PORTD = 0xAA;
    A2D_Setup();
    while(1){
        PORTD = A2D_getHighByteResult();
    }
}
void A2D_ReadA3ToPortD(){
    char a2dAnswer;
    TRISD = 0x00;
    PORTD = 0xAA; //Check lights
    // {00,0011,11} 0,Channel3,GO, ADON
    ADCON0 = 0x0D;
    // {0000,1011}=>0, Vref(+,-)=(AVdd, AVss), AN3-AN0 is analog
    ADCON1 = 0x0B;//TA
    //{0001,1110}=> LeftJustify,0, 6TAD, FOSC/64
    ADCON2 = 0x1E;//TA
    TRISAbits.RA3 =1; //Set the anlog channel to input
    while(1){
        a2dAnswer = A2D_getHighByteResult();
        PORTD = ADRESH;
    }
}

unsigned int A2D_readResult(){
    unsigned int result;
    char *pChar;
    pChar = (char *)&result;
    pChar[0] = ADRESL;
    pChar[1] = ADRESH;
    pChar[2] = 0;
    pChar[3] = 0;
    return result;
}
int binary_decimal(int n) /* Function to convert binary to decimal.*/
{
    int decimal=0, i=0, rem;
    while (n!=0) {
        rem = n%10; n/=10;
        decimal += rem*pow(2,i);
        ++i;
    }
    return decimal; }
//- See more at: http://www.programiz.com/c-programming/examples/binary-decimal-convert#sthash.1yr5g29e.dpuf
void A2D_Setup(){
    // Wrong Statements
//    ADCON0bits.ADON = 1; //Turn on ADC module
//    ADCON0bits.CHS = 0;
//    ADCON1bits.PCFG = 0b0000; // All Channels analog
//    ADCON1bits.VCFG = 0b00; // Vref(+,-)=(AVdd, AVss)
//    ADCON2bits.ADFM = 0; //Left Justified Result (ADRESHL[0:3]=0)
//    ADCON2bits.ACQT = 0b001; //2 Tad = A/D Acquisition Time Select Bits
//    ADCON2bits.ADCS = 0b001; // FOSC/8
            //0b111; // Frc (Clock derived from A/D RC oscillator)
//    ADCON0 = 0x0D;
//    ADCON1 = 0x0B;//TA    //0x0A;//JH
//    ADCON2 = 0x1E;//TA     //;0x14;//JH
//
//    TRISAbits.RA0 =1; //Set the anlog channel to input
//    TRISAbits.RA1 =1;
//    TRISAbits.RA2 =1;
//    TRISAbits.RA3 =1;
// {00,0011,11} 0,Channel3,GO, ADON

//ADCON0//0b00aaaabc (aaaa) Channel, b (Go/Done_L), c (ADON)
    ADCON0 = 0b00000101; //Channel(AN1),ADenabled
    //ADCON0 = 0x0D;//Channel(AN3), ADenabled
//ADCON1//0b00aabbbb (aa)Vref (bbbb) Port Config
    ADCON1 = 0x1B;//{0000,1011}=>0, Vref(+,-)=(AVdd, AVss), AN3-AN0 is analog
    //{0001,1110}=> LeftJustify,0, 6TAD, FOSC/64
    ADCON2 = 0x1E;//TA
    TRISAbits.RA3 =1; //Set the anlog channel to input


}
void A2D_start(){
    ADCON0bits.GO_NOT_DONE = 1;
}
int A2D_getBothByteResult(){
    A2D_start();
    int clockCycles=0;
    while(!ADCON0bits.GO_NOT_DONE){
        clockCycles++;
    }
    unsigned int result = A2D_readResult();
    result = result -15424;
    sprintf(A2D_ResultText, "%u", result);
    return result;
}
int A2D_getHighByteResult(){
    A2D_start();
    int clockCycles=0;
    while(ADCON0bits.GO_NOT_DONE){
        clockCycles++;
    }
    sprintf(A2D_ResultText, "%u", ADRESH);
    return ADRESH;
}

double     A2D_getTemperatureFromVolt(double voltage){
    //Assume 10mV/C or 19.5mV/C
    // VDD = 2.3-5.5 or 3.1 to 5.5
    // Vout = T_c*T+V_0C
    // Ta and V0C are defined in the data sheet. Units are 'mV/degreeCelsius' and 'mV'
    double THERMISTER_Ta =9.93;
    double THERMISTER_V0C =529;//550;// 510.0;
    //Ta=9.93, V0C 529 Good for 50-80
    //Ta=9.95, V0C 550 Good for 50-80
    double mv = voltage*1000;
    double T;
    T = (mv- THERMISTER_V0C)/THERMISTER_Ta; //(mv-500)/10;
    return T;
}
double  A2D_getTempFromA2D(){
    double V = A2D_getVoltFrom8bitA2D();
    //double V = A2D_getVoltFrom10bitA2D(ADRESH);
    return A2D_getTemperatureFromVolt(V);
}
double  A2D_getVoltFrom8bitA2D(){
    //double preGainOffset;
    double gain = 1;//2.9;//2.75;
    double Vrange = 1.61; //5
    return ADRESH/255.0*Vrange/gain;
}
double  A2D_getVoltFrom10bitA2D(){
    int fullA2D = A2D_readResult();
    //double preGainOffset;
    double gain = 1;//2.9;//2.75;
    double Vrange = 1.61; //5
    double A2Dmax = 255*4;
    return fullA2D/A2Dmax*Vrange/gain;
}
//        A2D functions         >>
//<<        Delay functions
void delayBetweenMeasure(){
    __delay_ms(50);
    __delay_ms(50);
}
void delay500ms(){
    int i = 0;
    for(i = 0;i<5;i++){
            __delay_ms(50);
            __delay_ms(50);
    }
    return;
}
void delayXSec(int loopNum)
{
    int i = 0;
    for(i = 0;i<loopNum*10;i++){
            __delay_ms(50);
            __delay_ms(50);
    }
    return;
}
//        Delay functions     >>