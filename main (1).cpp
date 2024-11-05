#include "mbed.h"
#include "rtos.h"
#include "EthernetInterface.h"
#include "ExperimentServer.h"
#include "QEI.h"
#include "BezierCurve.h"
#include "MotorShield.h" 
#include "HardwareSetup.h"

#define BEZIER_ORDER_FOOT_1    7
#define BEZIER_ORDER_FOOT_2    7
#define NUM_INPUTS (14 + 2*(BEZIER_ORDER_FOOT_1+1) + 2*(BEZIER_ORDER_FOOT_2+1))
#define NUM_OUTPUTS 38

#define PULSE_TO_RAD (2.0f*3.14159f / 1200.0f)

// Initializations
Serial pc(USBTX, USBRX);    // USB Serial Terminal
ExperimentServer server;    // Object that lets us communicate with MATLAB
Timer t;                    // Timer to measure elapsed time of experiment

QEI encoderA(PE_9,PE_11, NC, 1200, QEI::X4_ENCODING);  // MOTOR A ENCODER (no index, 1200 counts/rev, Quadrature encoding)
QEI encoderB(PA_5, PB_3, NC, 1200, QEI::X4_ENCODING);  // MOTOR B ENCODER (no index, 1200 counts/rev, Quadrature encoding)
QEI encoderC(PC_6, PC_7, NC, 1200, QEI::X4_ENCODING);  // MOTOR C ENCODER (no index, 1200 counts/rev, Quadrature encoding)
QEI encoderD(PD_12, PD_13, NC, 1200, QEI::X4_ENCODING);// MOTOR D ENCODER (no index, 1200 counts/rev, Quadrature encoding)

MotorShield motorShield(24000); //initialize the motor shield with a period of 24000 ticks or ~10kHZ
Ticker currentLoop;

// Variables for q1
float current1;
float current_des1 = 0;
float prev_current_des1 = 0;
float current_int1 = 0;
float angle1;
float velocity1;
float duty_cycle1;
float angle1_init;

// Variables for q2
float current2;
float current_des2 = 0;
float prev_current_des2 = 0;
float current_int2 = 0;
float angle2;
float velocity2;
float duty_cycle2;
float angle2_init;

// Variables for q3
float current3;
float current_des3 = 0;
float prev_current_des3 = 0;
float current_int3 = 0;
float angle3;
float velocity3;
float duty_cycle3;
float angle3_init;

// Variables for q4
float current4;
float current_des4 = 0;
float prev_current_des4 = 0;
float current_int4 = 0;
float angle4;
float velocity4;
float duty_cycle4;
float angle4_init;

// Fixed kinematic parameters
const float l_OA=.011; 
const float l_OB=.042; 
const float l_AC=.096; 
const float l_DE=.091;

// Timing parameters
float current_control_period_us = 200.0f;     // 5kHz current control loop
float impedance_control_period_us = 1000.0f;  // 1kHz impedance control loop
float start_period, traj_period, end_period;

// Control parameters
float current_Kp = 4.0f;         
float current_Ki = 0.4f;           
float current_int_max = 3.0f;       
float duty_max;      
float K_xx;
float K_yy;
float K_xy;
float D_xx;
float D_xy;
float D_yy;

// Model parameters
float supply_voltage = 12;     // motor supply voltage
float R = 2.0f;                // motor resistance
float k_t = 0.18f;             // motor torque constant
float nu = 0.0005;             // motor viscous friction

// Current control interrupt function
void CurrentLoop() {
    // This loop sets the motor voltage commands using PI current controllers with feedforward terms.
    
    //use the motor shield as follows:
    //motorShield.motorAWrite(DUTY CYCLE, DIRECTION), DIRECTION = 0 is forward, DIRECTION =1 is backwards.
        
    current1 = -(((float(motorShield.readCurrentA())/65536.0f)*30.0f)-15.0f);           // measure current
    velocity1 = encoderA.getVelocity() * PULSE_TO_RAD;                                  // measure velocity        
    float err_c1 = current_des1 - current1;                                             // current errror
    current_int1 += err_c1;                                                             // integrate error
    current_int1 = fmaxf( fminf(current_int1, current_int_max), -current_int_max);      // anti-windup
    float ff1 = R*current_des1 + k_t*velocity1;                                         // feedforward terms
    duty_cycle1 = (ff1 + current_Kp*err_c1 + current_Ki*current_int1)/supply_voltage;   // PI current controller
    
    float absDuty1 = abs(duty_cycle1);
    if (absDuty1 > duty_max) {
        duty_cycle1 *= duty_max / absDuty1;
        absDuty1 = duty_max;
    }    
    if (duty_cycle1 < 0) { // backwards
        motorShield.motorAWrite(absDuty1, 1);
    } else { // forwards
        motorShield.motorAWrite(absDuty1, 0);
    }             
    prev_current_des1 = current_des1; 
    
    current2 = -(((float(motorShield.readCurrentB())/65536.0f)*30.0f)-15.0f);           // measure current
    velocity2 = encoderB.getVelocity() * PULSE_TO_RAD;                                  // measure velocity  
    float err_c2 = current_des2 - current2;                                             // current error
    current_int2 += err_c2;                                                             // integrate error
    current_int2 = fmaxf( fminf(current_int2, current_int_max), -current_int_max);      // anti-windup   
    float ff2 = R*current_des2 + k_t*velocity2;                                         // feedforward terms
    duty_cycle2 = (ff2 + current_Kp*err_c2 + current_Ki*current_int2)/supply_voltage;   // PI current controller
    
    float absDuty2 = abs(duty_cycle2);
    if (absDuty2 > duty_max) {
        duty_cycle2 *= duty_max / absDuty2;
        absDuty2 = duty_max;
    }    
    if (duty_cycle2 < 0) { // backwards
        motorShield.motorBWrite(absDuty2, 1);
    } else { // forwards
        motorShield.motorBWrite(absDuty2, 0);
    }             
    prev_current_des2 = current_des2; 
    
    current3 = -(((float(motorShield.readCurrentA())/65536.0f)*30.0f)-15.0f);           // measure current
    velocity3 = encoderC.getVelocity() * PULSE_TO_RAD;                                  // measure velocity        
    float err_c3 = current_des1 - current1;                                             // current errror
    current_int3 += err_c3;                                                             // integrate error
    current_int3 = fmaxf( fminf(current_int3, current_int_max), -current_int_max);      // anti-windup
    float ff3 = R*current_des3 + k_t*velocity3;                                         // feedforward terms
    duty_cycle3 = (ff3 + current_Kp*err_c3 + current_Ki*current_int3)/supply_voltage;   // PI current controller
    
    float absDuty3 = abs(duty_cycle1);
    if (absDuty3 > duty_max) {
        duty_cycle3 *= duty_max / absDuty3;
        absDuty3 = duty_max;
    }    
    if (duty_cycle3 < 0) { // backwards
        motorShield.motorCWrite(absDuty3, 0);
    } else { // forwards
        motorShield.motorCWrite(absDuty3, 1);
    }             
    prev_current_des3 = current_des3; 
    
    current4 = -(((float(motorShield.readCurrentB())/65536.0f)*30.0f)-15.0f);           // measure current
    velocity4 = encoderD.getVelocity() * PULSE_TO_RAD;                                  // measure velocity  
    float err_c4 = current_des4 - current4;                                             // current error
    current_int4 += err_c4;                                                             // integrate error
    current_int4 = fmaxf( fminf(current_int4, current_int_max), -current_int_max);      // anti-windup   
    float ff4 = R*current_des4 + k_t*velocity4;                                         // feedforward terms
    duty_cycle4 = (ff4 + current_Kp*err_c4 + current_Ki*current_int4)/supply_voltage;   // PI current controller
    
    float absDuty4 = abs(duty_cycle4);
    if (absDuty4 > duty_max) {
        duty_cycle4 *= duty_max / absDuty4;
        absDuty4 = duty_max;
    }    
    if (duty_cycle4 < 0) { // backwards
        motorShield.motorDWrite(absDuty4, 0);
    } else { // forwards
        motorShield.motorDWrite(absDuty4, 1);
    }             
    prev_current_des4 = current_des4; 
}

// JUST TO CHECK
// Directly set motor speeds without current control loop
void DirectMotorControl() {
    // Set a constant duty cycle for all motors to make them run at the same pace
    float constant_duty_cycle = 0.5f; // Example duty cycle value between 0 and 1

    // Set motor A (forward direction)
    motorShield.motorAWrite(constant_duty_cycle, 0); // 0 = forward direction

    // Set motor B (forward direction)
    motorShield.motorBWrite(constant_duty_cycle, 0); // 0 = forward direction

    // Set motor C (forward direction)
    motorShield.motorCWrite(constant_duty_cycle, 0); // 0 = forward direction

    // Set motor D (forward direction)
    motorShield.motorDWrite(constant_duty_cycle, 0); // 0 = forward direction
}

int main(void) {
    // Object for 7th order Cartesian foot trajectory
    BezierCurve rDesFoot_bez_1(2,BEZIER_ORDER_FOOT_1);
    BezierCurve rDesFoot_bez_2(2,BEZIER_ORDER_FOOT_2);
    
    // Link the terminal with our server and start it up
    server.attachTerminal(pc);
    server.init();
    
    // Continually get input from MATLAB and run experiments
    float input_params[NUM_INPUTS];
    pc.printf("%f",input_params[0]);
    
    while(1) {
        // If there are new inputs, this code will run
        if (server.getParams(input_params,NUM_INPUTS)) {          
            // Get inputs from MATLAB          
            start_period                = input_params[0];    // First buffer time, before trajectory
            traj_period                 = input_params[1];    // Trajectory time/length
            end_period                  = input_params[2];    // Second buffer time, after trajectory
    
            angle1_init                 = input_params[3];    // Initial angle for q1 (rad)
            angle2_init                 = input_params[4];    // Initial angle for q2 (rad)
            angle3_init                 = input_params[12];    // Initial angle for q3 (rad)
            angle4_init                 = input_params[13];    // Initial angle for q4 (rad)

            K_xx                        = input_params[5];    // Foot stiffness N/m
            K_yy                        = input_params[6];    // Foot stiffness N/m
            K_xy                        = input_params[7];    // Foot stiffness N/m
            D_xx                        = input_params[8];    // Foot damping N/(m/s)
            D_yy                        = input_params[9];    // Foot damping N/(m/s)
            D_xy                        = input_params[10];   // Foot damping N/(m/s)
            duty_max                    = input_params[11];   // Maximum duty factor
          
            // Get foot trajectory points
            float foot_pts1[2*(BEZIER_ORDER_FOOT_1+1)];
            for(int i = 0; i<2*(BEZIER_ORDER_FOOT_1+1);i++) {
              foot_pts1[i] = input_params[14+i];    
            }
            rDesFoot_bez_1.setPoints(foot_pts1);

            float foot_pts2[2*(BEZIER_ORDER_FOOT_2+1)];
            for(int i = 0; i<2*(BEZIER_ORDER_FOOT_2+1);i++) {
              foot_pts2[i] = input_params[14+i];    
            }
            rDesFoot_bez_2.setPoints(foot_pts2);
            
            // Attach current loop interrupt
            currentLoop.attach_us(CurrentLoop,current_control_period_us);

            // JUST TO CHECK
            // Start motors at a constant duty cycle
            // DirectMotorControl();
                        
            // Setup experiment
            t.reset();
            t.start();
            encoderA.reset();
            encoderB.reset();
            encoderC.reset();
            encoderD.reset();

            motorShield.motorAWrite(0, 0); //turn motor A off
            motorShield.motorBWrite(0, 0); //turn motor B off
            motorShield.motorCWrite(0, 0); //turn motor C off
            motorShield.motorDWrite(0, 0); //turn motor D off
                         
            // Run experiment
            while( t.read() < start_period + traj_period + end_period) { 
                DirectMotorControl();

                // Read encoders to get motor states
                angle1 = encoderA.getPulses() *PULSE_TO_RAD + angle1_init;       
                velocity1 = encoderA.getVelocity() * PULSE_TO_RAD;
                 
                angle2 = encoderB.getPulses() * PULSE_TO_RAD + angle2_init;       
                velocity2 = encoderB.getVelocity() * PULSE_TO_RAD;   

                angle3 = encoderC.getPulses() *PULSE_TO_RAD + angle3_init;       
                velocity3 = encoderC.getVelocity() * PULSE_TO_RAD;
                 
                angle4 = encoderD.getPulses() * PULSE_TO_RAD + angle4_init;       
                velocity4 = encoderD.getVelocity() * PULSE_TO_RAD;           
                
                const float th1 = angle1;
                const float th2 = angle2;
                const float th3 = angle3;
                const float th4 = angle4;
                const float dth1= velocity1;
                const float dth2= velocity2;
                const float dth3= velocity3;
                const float dth4= velocity4;
 
                // Calculate the Jacobian
                float Jx_th1 = l_OB * cos(th1) + l_AC * cos(th1 + th2) + l_DE * cos(th1);
                float Jy_th1 = l_OB * sin(th1) + l_AC * sin(th1 + th2) + l_DE * sin(th1);
                float Jx_th2 = l_AC * cos(th1 + th2);
                float Jy_th2 = l_AC * sin(th1 + th2);
                                
                // Calculate the forward kinematics (position and velocity)
                float xFoot1 = l_OB * sin(th1) + l_AC * sin(th1 + th2) + l_DE * sin(th1);
                float yFoot1 = -l_OB * cos(th1) - l_AC * cos(th1 + th2) - l_DE * cos(th1);
                float dxFoot1 = l_OB * cos(th1) * dth1 + l_AC * cos(th1 + th2) * (dth1 + dth2) + l_DE * cos(th1) * dth1; 
                float dyFoot1 = l_OB * sin(th1) * dth1 + l_AC * sin(th1 + th2) * (dth1 + dth2) + l_DE * sin(th1) * dth1;  
                
                float xFoot2 = l_OB * sin(th3) + l_AC * sin(th3 + th4) + l_DE * sin(th3);
                float yFoot2 = -l_OB * cos(th3) - l_AC * cos(th3 + th4) - l_DE * cos(th3);
                float dxFoot2 = l_OB * cos(th3) * dth1 + l_AC * cos(th3 + th4) * (dth3 + dth4) + l_DE * cos(th3) * dth3; 
                float dyFoot2 = l_OB * sin(th3) * dth1 + l_AC * sin(th3 + th4) * (dth3 + dth4) + l_DE * sin(th3) * dth3; 
                    


                // Set gains based on buffer and traj times, then calculate desired x,y from Bezier trajectory at current time if necessary
                float teff  = 0;
                float vMult = 0;
                if (t < start_period) {
                    if (K_xx > 0 || K_yy > 0) {
                        K_xx = 1; // for joint space control, set this to 1; for Cartesian space control, set this to 50
                        K_yy = 1; // for joint space control, set this to 1; for Cartesian space control, set this to 50
                        D_xx = 0.1;  // for joint space control, set this to 0.1; for Cartesian space control, set this to 2
                        D_yy = 0.1;  // for joint space control, set this to 0.1; for Cartesian space control, set this to 2
                        K_xy = 0;
                        D_xy = 0;
                    }
                    teff = 0;
                }
                else if (t < start_period + traj_period) {
                    K_xx = input_params[5];  // Foot stiffness N/m
                    K_yy = input_params[6];  // Foot stiffness N/m
                    K_xy = input_params[7];  // Foot stiffness N/m
                    D_xx = input_params[8];  // Foot damping N/(m/s)
                    D_yy = input_params[9];  // Foot damping N/(m/s)
                    D_xy = input_params[10]; // Foot damping N/(m/s)
                    teff = (t-start_period);
                    vMult = 1;
                }
                else {
                    teff = traj_period;
                    vMult = 0;
                }
                
                // Get desired foot positions and velocities
                float rDesFoot1[2], vDesFoot1[2];
                rDesFoot_bez_1.evaluate(teff/traj_period,rDesFoot1);
                rDesFoot_bez_1.evaluateDerivative(teff/traj_period,vDesFoot1);
                vDesFoot1[0]/=traj_period;
                vDesFoot1[1]/=traj_period;
                vDesFoot1[0]*=vMult;
                vDesFoot1[1]*=vMult;
                
                // Calculate the inverse kinematics (joint positions and velocities) for desired joint angles              
                float xFoot1_inv = -rDesFoot1[0];
                float yFoot1_inv = rDesFoot1[1];                
                float l_OE_1 = sqrt( (pow(xFoot1_inv,2) + pow(yFoot1_inv,2)) );
                float alpha = abs(acos( (pow(l_OE_1,2) - pow(l_AC,2) - pow((l_OB+l_DE),2))/(-2.0f*l_AC*(l_OB+l_DE)) ));
                float th2_des_1 = -(3.14159f - alpha); 
                float th1_des_1 = -((3.14159f/2.0f) + atan2(yFoot1_inv,xFoot1_inv) - abs(asin( (l_AC/l_OE_1)*sin(alpha) )));

                float rDesFoot2[2], vDesFoot2[2];
                rDesFoot_bez_2.evaluate(teff/traj_period,rDesFoot2);
                rDesFoot_bez_2.evaluateDerivative(teff/traj_period,vDesFoot2);
                vDesFoot2[0]/=traj_period;
                vDesFoot2[1]/=traj_period;
                vDesFoot2[0]*=vMult;
                vDesFoot2[1]*=vMult;

                float xFoot2_inv = -rDesFoot2[0];
                float yFoot2_inv = rDesFoot2[1];                
                float l_OE_2 = sqrt( (pow(xFoot2_inv,2) + pow(yFoot2_inv,2)) );
                float th2_des_2 = -(3.14159f - alpha); 
                float th1_des_2 = -((3.14159f/2.0f) + atan2(yFoot2_inv,xFoot2_inv) - abs(asin( (l_AC/l_OE_2)*sin(alpha) )));
                
                float dd = (Jx_th1*Jy_th2 - Jx_th2*Jy_th1);
                float dth1_des = (1.0f/dd) * (  Jy_th2*vDesFoot1[0] - Jx_th2*vDesFoot1[1] );
                float dth2_des = (1.0f/dd) * ( -Jy_th1*vDesFoot1[0] + Jx_th1*vDesFoot1[1] );
                float dth3_des = (1.0f/dd) * (  Jy_th2*vDesFoot2[0] - Jx_th2*vDesFoot2[1] );
                float dth4_des = (1.0f/dd) * ( -Jy_th1*vDesFoot2[0] + Jx_th1*vDesFoot2[1] );
        
                // Calculate error variables
                float e_x_1 = rDesFoot1[0] - xFoot1;
                float e_y_1 = rDesFoot1[1] - yFoot1;
                float de_x_1 = vDesFoot1[0] - dxFoot1;
                float de_y_1 = vDesFoot1[1] - dyFoot1;
                float e_x_2 = rDesFoot2[0] - xFoot2;
                float e_y_2 = rDesFoot2[1] - yFoot2;
                float de_x_2 = vDesFoot2[0] - dxFoot2;
                float de_y_2 = vDesFoot2[1] - dyFoot2;
        
                // Calculate virtual force on foot
                float fx = 0;
                float fy = 0;
                                
                // Set desired currents             
                current_des1 = 0;          
                current_des2 = 0;   
        
                // Joint impedance
                // sub Kxx for K1, Dxx for D1, Kyy for K2, Dyy for D2
                // Note: Be careful with signs now that you have non-zero desired angles!
                // Your equations should be of the form i_d = K1*(q1_d - q1) + D1*(dq1_d - dq1)
//                current_des1 = 0;          
//                current_des2 = 0;                          
                           
                // Cartesian impedance  
                // Note: As with the joint space laws, be careful with signs!              
//                current_des1 = 0;          
//                current_des2 = 0;   
                
                
                // Form output to send to MATLAB     
                float output_data[NUM_OUTPUTS];
                output_data[0] = t.read();
                // motor 1 state
                output_data[1] = angle1;
                output_data[2] = velocity1;  
                output_data[3] = current1;
                output_data[4] = current_des1;
                output_data[5] = duty_cycle1;
                // motor 2 state
                output_data[6] = angle2;
                output_data[7] = velocity2;
                output_data[8] = current2;
                output_data[9] = current_des2;
                output_data[10]= duty_cycle2;
                // motor 3 state
                output_data[11] = angle2;
                output_data[12] = velocity2;
                output_data[13] = current2;
                output_data[14] = current_des2;
                output_data[15]= duty_cycle2;
                // motor 4 state
                output_data[16] = angle2;
                output_data[17] = velocity2;
                output_data[18] = current2;
                output_data[19] = current_des2;
                output_data[20]= duty_cycle2;
                // foot 1 state
                output_data[21] = xFoot1;
                output_data[22] = yFoot1;
                output_data[23] = dxFoot1;
                output_data[24] = dyFoot1;
                output_data[25] = rDesFoot1[0];
                output_data[26] = rDesFoot1[1];
                output_data[27] = vDesFoot1[0];
                output_data[28] = vDesFoot1[1];
                //foot 2 state
                output_data[29] = xFoot2;
                output_data[30] = yFoot2;
                output_data[32] = dxFoot2;
                output_data[33] = dyFoot2;
                output_data[34] = rDesFoot2[0];
                output_data[35] = rDesFoot2[1];
                output_data[36] = vDesFoot2[0];
                output_data[37] = vDesFoot2[1];
                // Send data to MATLAB
                server.sendData(output_data,NUM_OUTPUTS);
                
                // Send data to MATLAB
                server.sendData(output_data,NUM_OUTPUTS);

                wait_us(impedance_control_period_us);   
            }
            
            // Cleanup after experiment
            server.setExperimentComplete();
            currentLoop.detach();
            motorShield.motorAWrite(0, 0); //turn motor A off
            motorShield.motorBWrite(0, 0); //turn motor B off
            motorShield.motorCWrite(0, 0); //turn motor A off
            motorShield.motorDWrite(0, 0); //turn motor B off
        } // end if
    } // end while
} // end main

