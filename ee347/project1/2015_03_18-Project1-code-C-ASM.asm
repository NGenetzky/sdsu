;/*
; * File:   main.c
; * Author: Nathansen
; *
; * PORTA is used for Anlog Sensor
; *  A1 Thermister Reading
; *  A3 Vref+
; * PORTE is used for keypad.
; *  E0  COL0
; *  E1  COL1
; *  E2  COL2
; *  E3  COL3
; *  E4  ROW0
; *  E5  ROW1
; *  E6  ROW2
; *  E7  ROW3
; * PORTJ is used for the TEA input
; *  J0  Run Switch
; *  J1  Input Switch
; * PORTJ is used for LED.
; *  J3  Green LED
; *  J4  Amber LED
; *  J5  Red LED
; *  J6  Run LED
; *  J7  Heater LED
; * Created on January 24, 2015, 3:13 PM
; */
 LIST P=18f8722, F=INHX32, N=0, ST=OFF, R=DEC

#include <p18F8722.inc>
        ;CONFIG OSC=HS
        CONFIG OSC = HSPLL    ; Oscillator Selection bits (HS oscillator)
        CONFIG WDT=OFF
        CONFIG LVP = OFF   ;//Single-Supply ISCP Enable bit
DELAYTIME EQU 0x43
DELAYREG3 EQU 0x42
DELAYREG2 EQU 0x41
DELAYREG1 EQU 0x40
DELAYHTIME EQU 0x47
DELAYHREG3 EQU 0x46
DELAYHREG2 EQU 0x45
DELAYHREG1 EQU 0x44
keypad		equ	PORTE		; keypad port
keypd_dir	equ	TRISE		; keypad direction register

;This code is to convert hex to string on the keypad
;KEY00       equ '0'
;KEY01       equ '1'
;KEY02       equ '2'
;KEY03       equ '3'
KEY00       equ D'1'
KEY01       equ D'2'
KEY02       equ D'4'
KEY03       equ D'8'
KEY10       equ D'4'
KEY11       equ D'5'
KEY12       equ D'6'
KEY13       equ D'7'
KEY20       equ D'8'
KEY21       equ D'9'
KEY22       equ D'10'
KEY23       equ D'11'
KEY30       equ D'12'
KEY31       equ D'13'
KEY32       equ D'14'
KEY33       equ D'15'
;OUTPUT
#define TEA_TRIS 0x03;
#define PORT_GREENLED PORTJ,3 ; To Use: "BSF PORT_GREENLED" or "BCF PORT_GREENLED"
#define PORT_AMBERLED PORTJ,4
#define PORT_REDLED PORTJ,5
#define PORT_RUNLED PORTJ,6
#define PORT_HEATER PORTJ,7
;INPUT
#define PORT_SWITCH_RUN PORTJ,0
#define PORT_SWITCH_INPUT PORTJ,1

;Temperatures
;0.6612V=106/255=31.5C
#define TEMP23 116
#define TEMP55 165
#define TEMP78 205
#define TEMP80 207
#define TEMP85 213
#define TEMP92 223
#define TEMP95 232

;VARIABLE DECLARATIONS
    UDATA_ACS
;UDATA_ACS************************************************************* TEA
low_temp res 1    ;initialized
high_temp res 1   ;initialized
temperature res 1 ;initialized

STATE_FULLBOIL_H  equ TEMP95+4
STATE_FULLBOIL_L  equ TEMP95+2
STATE_SIMMER_H    equ TEMP95+1
STATE_SIMMER_L    equ TEMP95-1
STATE_SOUP_H      equ TEMP80+3
STATE_SOUP_L      equ TEMP80-2
;UDATA_ACS************************************************************* TEA
;UDATA_ACS **************************************************************KEYPAD
R1      res 1
R2      res 1               
R3      res 1
dtWREG  res 1               ; so delay subR can preserve WREG
alupW   res 1               ; make all_up save WREG also
alupLV  res 1               ; do "allup multiple times" w/debounce in between?
;UDATA_ACS **************************************************************KEYPAD
		org	0x00
		goto	start
		org	0x08
		retfie
		org	0x18
		retfie
start
        clrf    TRISD       
        clrf    PORTD
        MOVLW TEA_TRIS
        MOVWF TRISJ
        MOVLW D'105'
        MOVWF low_temp
        MOVLW D'120'
        MOVWF high_temp
        ;call    A2D_AN1_PosVRef   ;calls setup_a2d subroutine
;Main loop for the program
disabled
    BCF PORT_RUNLED
    BCF PORT_GREENLED   ; Turn on Green Light
    BCF PORT_REDLED
    BCF PORT_AMBERLED
    BCF PORT_HEATER
loop
        ;call TwoFtyMS_DELAY
        ;call GET_A2D
        ;MOVWF PORTD; Read the high part of A2D
        ;BRA loop

        ;BTFSS PORT_SWITCH_RUN
            ;BRA disabled ;Infinite loop if RUN switch is off
        BTFSC PORT_SWITCH_INPUT
            BRA USER_INPUT ;go to UserInput if INPUT switch is on
        BRA Feedbackloop ;go to Feedbackloop if INPUT switch is off
		
        BRA	loop ;Never Reach
		goto	$ ;Loop forever

;**************************************************************  USER INPUT LOOP
; USER_INPUT will occur if the INPUT switch is on
; Will wait for key press. First 3 buttons set the state of TEA
; Will return to the main loop after
USER_INPUT
        call delay10ms
        BTFSS PORT_SWITCH_INPUT ; If switch is OFF then go back to loop
            BRA loop ; debounce
        call get_key
        MOVWF PORTD
        btfsc   WREG,0
            bra SET_STATE_SOUP
        btfsc   WREG,1    ;
            bra SET_STATE_SIMMER
        btfsc   WREG,2
            bra SET_STATE_FULLBOIL
        btfsc   WREG,3
            bra loop ; Button 3 has no function
        bra loop ; Button 4-15 have no function
;sets to FULLBOIL state (100 degrees celcius)
SET_STATE_FULLBOIL
;     BSF    PORTD,5
;     BCF    PORTD,6
;     BCF    PORTD,7
     MOVLW STATE_FULLBOIL_H
     MOVWF high_temp
     MOVLW STATE_FULLBOIL_L
     MOVWF low_temp
     BRA loop ; Return to main loop

;sets to SIMMER state (90 degrees celcius)
SET_STATE_SIMMER
;     BCF    PORTD,5
;     BSF    PORTD,6
;     BCF    PORTD,7
     MOVLW STATE_SIMMER_H
     MOVWF high_temp
     MOVLW STATE_SIMMER_L
     MOVWF low_temp
     BRA loop ; Return to main loop

;State SOUP 80degrees celcius
SET_STATE_SOUP ;sets the soup state with the keypad
;     BCF    PORTD,5
;     BCF    PORTD,6
;     BSF    PORTD,7

     MOVLW STATE_SOUP_H
     MOVWF high_temp
     MOVLW STATE_SOUP_L
     MOVWF low_temp
     BRA loop ; Return to main loop
;**************************************************************  USER INPUT LOOP
;****************************************************************  FEEDBACK LOOP
; Feedbackloop will occur if the INPUT switch is off
; Will wait for A2D to read. Will set LEDs and Heater
; ...based on current reading vs desired reading
; Will return to the main loop after
Feedbackloop
    call delay10ms
    BTFSC PORT_SWITCH_INPUT
        BRA loop ; If switch is ON then go back to loop
    BSF PORT_RUNLED
    call GET_A2D    ; Get value of the A2D
    MOVWF PORTD
    MOVWF temperature
    MOVF high_temp, w
    SUBWF temperature,w ;Subtract temperature-WREG
    BNN temp_is_high ; if (temperature-high_temp >0)
    MOVF low_temp, w
    SUBWF temperature,w
    BNN temp_is_right ; if (temp-lowTemp >0)
    BRA temp_is_low ; if (temp-lowTemp <0)

;If temperature too high turn off heater to lower temp
temp_is_high
    BSF PORT_REDLED ; Turn on Red Light
    BCF PORT_GREENLED
    BCF PORT_AMBERLED
    BCF PORT_HEATER
    BRA loop ; Return to main loop

;If temperture is right
temp_is_right
    BSF PORT_GREENLED   ; Turn on Green Light
    BCF PORT_REDLED
    BCF PORT_AMBERLED
    BCF PORT_HEATER
    BRA loop ; Return to main loop

;If temperature is too low turn on heater to raise temp
temp_is_low
    BSF PORT_AMBERLED   ; Turn on Amber Light
    BCF PORT_GREENLED
    BCF PORT_REDLED
    BSF PORT_HEATER     ; Turn on Heater
    BRA loop ; Return to main loop
;****************************************************************  FEEDBACK LOOP
;************************************************************************  A2D
;The set up for A2D
setup_a2d
 CLRF TRISD
    BSF TRISA, 3
    MOVLW B'00001101' ;AN1
    MOVWF ADCON0
    MOVLW B'00001010'  ;Left Justified
    MOVWF ADCON1
    MOVLW B'00010100';6TAD, FOSC/4
    MOVWF ADCON2
    return
;end of set up for A2D
;The set up for A2D
A2D_AN1_PosVRef
    BSF TRISA, 1
    BSF TRISA, 3 
    MOVLW B'00000101' ;Channel(AN1),ADenabled
    MOVWF ADCON0
    MOVLW 0x1B ;MOVLW B'00011011'  ;{0000,1011}=>0, Vref(+,-)=(AVdd, AVss), AN3-AN0 is analog
    MOVWF ADCON1
    MOVLW B'00011110' ;{0001,1110}=> LeftJustify,0, 6TAD, FOSC/64
    MOVWF ADCON2
    return
;end of set up for A2D

;get reading from a2d
GET_A2D
    call A2D_AN1_PosVRef
    BSF ADCON0, 1;GO
BACK ;loop to check if done
    BTFSC ADCON0, 1;DONE
    BRA BACK ; End loop to check if done
    ;MOVFF ADRESH, PORTD ;Moves the answer. Want in WREG.
    MOVF ADRESH,W
    return
;************************************************************************  A2D
;********************************************************************   KEYPAD
;*********************************
;
;       This code is provded to EE 347 students in the
;       Spring 2015 semester.  It will be considerd an act
;       of academic dishonesty to either pass this code alont
;       to students taking this class at a later date or
;       to use this code at a later date w/out specific
;       permision from the course instructor.
;
;       make port easy to change for demo
;       (and in case a damaged board is encountered)
;
;       Students should try Port H  first, and look at the
;       output using the digital capabilities of the scopes
;       in DEH 232.  A breakpoint should be set at the line:
;       "bcf	keypad,4"   below
;       which will show that PORTH  is not a good idea, since
;       clearing a singe bit clears (at least) the entire
;       high nibble.  (I still need to google to see why this
;       is the case, but this happened on several boards)
;
;       I then tried several ports and this code ended up
;       working if PORTC was used (change both PORT_ and TRIS_
;       below
;
;*********************************************
;
; SubR all_up   Check that no key is being pushed
;
; Quick and dirty for debugging race condition.
; won't work if you hold one above down before
; releasing one below it...
;
; If this stays in it should be improved.  Turned out to
; not be needed, so it's rough but not called
;
;**************************
all_up
        movwf   alupW
        movlw	0x0F		; configure the upper four pins of keypad
		movwf	keypd_dir 	; port for output, others for input
        movlw   0x0F
 		andwf	keypad,F    ; drive 0 on all rows
        movlw   2           ; debug, do this twice?
        movwf   alupLV
c0dwn   btfss   keypad,0
        bra c0dwn           ;just scans one at a time, doesn't recheck @ end
c1dwn   btfss   keypad,1    ; (won't work if you hold one down before releasing the previous)
        bra c1dwn
c2dwn   btfss   keypad,2
        bra c2dwn
c3dwn   btfss   keypad,3
        bra c3dwn
        call delay10ms
        decfsz  alupLV
        bra c0dwn           ; debounce and do 2x
        movlw	0xF0        ; Drive ones on all rows again (what get_key needs/does
	iorwf	keypad,F
        nop
        movf alupW,W
        return
; ***********************************************************************
; scan keypad, if hit, wait to debounce, and return
; the character to the calling routine.
; Called by driver routine for now,
; ***********************************************************************
get_key	movlw	0x0F		; configure the upper four pins of keypad
		movwf	keypd_dir 	; port for output, others for input
        movlw	0xF0		;
		iorwf	keypad,F
scan_r0	movlw	0xFF		; prepare to scan row 0 (driven low by RH4)
		iorwf	keypad,F 	; 	"
		bcf	keypad,4        ; students should set breakpoint here
scan_k0	btfss	keypad,0 	; check key 0
		goto 	db_key0
scan_k1	btfss	keypad,1 	; check key 1
		goto	db_key1
scan_k2	btfss	keypad,2 	; check key 2
		goto	db_key2
scan_k3	btfss	keypad,3 	; check key 3
		goto	db_key3
scan_r1	movlw	0xFF		; prepare to scan row 1 (driven low by RH5)
		iorwf	keypad,F
		bcf	keypad,5
scan_k4	btfss	keypad,0 	; check key 4
		goto	db_key4
scan_k5	btfss	keypad,1 	; check key 5
		goto	db_key5
scan_k6	btfss	keypad,2 	; check key 6
		goto	db_key6
scan_k7	btfss	keypad,3 	; check key 7
		goto	db_key7
scan_r2	movlw	0xFF		; prepare to scan row 2 (driven low by RH6)
		iorwf	keypad,F
		bcf	keypad,6
scan_k8	btfss	keypad,0 	; check key 8
		goto	db_key8
scan_k9	btfss	keypad,1 	; check key 9
		goto	db_key9
scan_kA	btfss	keypad,2 	; check key A
		goto	db_keyA
scan_kB	btfss	keypad,3 	; check key B
		goto	db_keyB
scan_r3	movlw	0xFF		; prepare to scan row 3 (driven low by RH7)
		iorwf	keypad,F	; 	"
		bcf	keypad,7        ;
scan_kC	btfss	keypad,0	; check key C
		goto	db_keyC
scan_kD	btfss	keypad,1	; check key D
		goto	db_keyD
scan_kE	btfss	keypad,2	; check key E
		goto	db_keyE
scan_kF	btfss	keypad,3	; check key F
		goto	db_keyF
		goto	scan_r0
db_key0	call	delay10ms	; wait for 10 ms to debounce
		btfsc	keypad,0
		goto	scan_k1     ; key 0 not pressed, check key 1
		movlw	KEY00  ;0x30
		return
db_key1	call	delay10ms
		btfsc	keypad,1
		goto	scan_k2
		movlw	KEY01 ;0x31
		return
db_key2	call	delay10ms
		btfsc	keypad,2
		goto	scan_k3
		movlw	KEY02 ;0x32
		return
db_key3	call	delay10ms
		btfsc	keypad,3
		goto	scan_r1
		movlw	KEY03 ;0x33
		return
db_key4	call	delay10ms
		btfsc	keypad,0
		goto	scan_k5
		movlw	KEY10 ;0x34
		return
db_key5	call	delay10ms
		btfsc	keypad,1
		goto	scan_k6
		movlw	KEY11 ;0x35
		return
db_key6	call	delay10ms
		btfsc	keypad,2
		goto	scan_k7
		movlw	KEY12 ;0x36
		return
db_key7	call	delay10ms
		btfsc	keypad,3
		goto	scan_r2
		movlw	KEY13 ;0x37
		return
db_key8	call	delay10ms
		btfsc	keypad,0
		goto	scan_k9
		movlw	KEY20 ;0x38
		return
db_key9	call	delay10ms
		btfsc	keypad,1
		goto	scan_kA
		movlw	KEY21 ;0x39
		return
db_keyA	call	delay10ms
		btfsc	keypad,2
		goto	scan_kB
		movlw	KEY22 ;KEY0x41
		return
db_keyB	call	delay10ms
		btfsc	keypad,3
		goto	scan_r3
		movlw	KEY23 ;0x42
		return
db_keyC	call	delay10ms
		btfsc	keypad,0
		goto	scan_kD
		movlw	KEY30 ;0x43
		return
db_keyD	call	delay10ms
		btfsc	keypad,1
		goto	scan_kE
		movlw	KEY31 ;0x44
		return
db_keyE	call	delay10ms
		btfsc	keypad,2
		goto	scan_kF
		movlw	KEY32 ;0x45
		return
db_keyF	call	delay10ms
		btfsc	keypad,3
		goto	scan_r0
		movlw	KEY33 ;0x46
		return
; KEYPAD_GetKey() DONE

;Delay subroutines for both A2D and KEYPAD CODE

;**********Delay subroutine for KEYPAD
; Use timers next week, for now use
; delay loop to wait 10ms for debounce
;*************************
delay10ms
        nop
        movwf dtWREG  ; save WREG, let it trash "STATUS" for now
        movlw D'1'    ; debug loop, needs to be one for sane value
        movwf R1
D1      movlw D'25'   ; should be d25 for 10 ms
        movwf R2
D2      movlw D'250'
        movwf  R3
D3      nop      ; Waste clock pulse
        nop
        decf R3, F
        bnz  D3
        decf  R2, F
        bnz  D2
        decf R1, F
        bnz D1
        movf dtWREG, W ; restore WREG (which may have key value)
        return


 ;This is the A2D Delay code
TwoFtyMS_DELAY
Q0
	MOVLW 0x0D	;Loads WREG with #
	MOVWF DELAYREG1	; Moves WREG To Delay REG 1
Q1
	MOVLW 0xFF	;Loads WREG with #
	MOVWF DELAYREG2	; Moves WREG To Delay REG 2
Q2
	MOVLW 0x26	;Loads WREG with #
	MOVWF DELAYREG3	;Moves WREG To Delay REG 3
Q3
	NOP
	NOP
	DECF DELAYREG3, F	;Decrements DELAY Reg 3
	BNZ Q3				;If Zero Flag bit is 0, branch to Q3
	DECF DELAYREG2, F	;Decrements DELAY Reg 3
	BNZ Q2				;If Zero Flag bit is 0, branch to Q2
	DECF DELAYREG1, F	;Decrements DELAY Reg 3
	BNZ Q1				;If Zero Flag bit is 0, branch to Q1
	DECF DELAYTIME, F
	BNZ Q0
    RETURN

    END