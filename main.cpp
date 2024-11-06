#include "mbed.h"
#include "rtos.h"
#include "EthernetInterface.h"
#include "ExperimentServer.h"
#include "QEI.h"
#include "BezierCurve.h"
#include "MotorShield.h" 
#include "HardwareSetup.h"
#include "Matrix.h"
#include "MatrixMath.h"

#define BEZIER_ORDER_FOOT    7
#define NUM_INPUTS (12 + 2*(BEZIER_ORDER_FOOT+1))
#define NUM_OUTPUTS 19

#define PULSE_TO_RAD (2.0f*3.14159f / 1200.0f)

// Initializations
Serial pc(USBTX, USBRX);    // USB Serial Terminal
ExperimentServer server;    // Object that lets us communicate with MATLAB
Timer t;                    // Timer to measure elapsed time of experiment

QEI encoderA(PE_9,PE_11, NC, 1200, QEI::X4_ENCODING);  // MOTOR A ENCODER (no index, 1200 counts/rev, Quadrature encoding)
QEI encoderB(PA_5, PB_3, NC, 1200, QEI::X4_ENCODING);  // MOTOR B ENCODER (no index, 1200 counts/rev, Quadrature encoding)
QEI encoderC(PC_6, PC_7, NC, 1200, QEI::X4_ENCODING);  // MOTOR C ENCODER (no index, 1200 counts/rev, Quadrature encoding)
QEI encoderD(PD_12, PD_13, NC, 1200, QEI::X4_ENCODING);// MOTOR D ENCODER (no index, 1200 counts/rev, Quadrature encoding)

MotorShield motorShield(24000); //initialize the motor shield with a period of 12000 ticks or ~20kHZ
Ticker currentLoop;

Matrix MassMatrix_1(2,2);
Matrix Jacobian_1(2,2);
Matrix JacobianT_1(2,2);
Matrix InverseMassMatrix_1(2,2);
Matrix temp_product_1(2,2);
Matrix Lambda_1(2,2);

Matrix MassMatrix_2(2,2);
Matrix Jacobian_2(2,2);
Matrix JacobianT_2(2,2);
Matrix InverseMassMatrix_2(2,2);
Matrix temp_product_2(2,2);
Matrix Lambda_2(2,2);

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
const float m1 =.0393 + .2;
const float m2 =.0368; 
const float m3 = .00783;
const float m4 = .0155;
const float I1 = 0.0000251;  //25.1 * 10^-6;
const float I2 = 0.0000535;  //53.5 * 10^-6;
const float I3 = 0.00000925; //9.25 * 10^-6;
const float I4 = 0.0000222;  //22.176 * 10^-6;
const float l_O_m1=0.032;
const float l_B_m2=0.0344; 
const float l_A_m3=0.0622;
const float l_C_m4=0.0610;
const float N = 18.75;
const float Ir = 0.0035/pow(N,2);

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
void CurrentLoop()
{
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
    
    current2     = -(((float(motorShield.readCurrentB())/65536.0f)*30.0f)-15.0f);       // measure current
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

    current3 = -(((float(motorShield.readCurrentC())/65536.0f)*30.0f)-15.0f);           // measure current
    velocity3 = encoderC.getVelocity() * PULSE_TO_RAD;                                  // measure velocity        
    float err_c3 = current_des3 - current3;                                             // current errror
    current_int3 += err_c3;                                                             // integrate error
    current_int3 = fmaxf( fminf(current_int3, current_int_max), -current_int_max);      // anti-windup
    float ff3 = R*current_des3 + k_t*velocity3;                                         // feedforward terms
    duty_cycle3 = (ff3 + current_Kp*err_c3 + current_Ki*current_int3)/supply_voltage;   // PI current controller
    
    float absDuty3 = abs(duty_cycle3);
    if (absDuty3 > duty_max) {
        duty_cycle3 *= duty_max / absDuty3;
        absDuty3 = duty_max;
    }    
    if (duty_cycle3 < 0) { // backwards
        motorShield.motorCWrite(absDuty3, 1);
        //motorShield.motorCWrite(absDuty1, 0);
    } else { // forwards
        motorShield.motorCWrite(absDuty3, 0);
        //motorShield.motorCWrite(absDuty1, 1);
    }             
    prev_current_des3 = current_des3; 
    

    current4     = -(((float(motorShield.readCurrentD())/65536.0f)*30.0f)-15.0f);       // measure current
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
        motorShield.motorDWrite(absDuty4, 1);
        //motorShield.motorDWrite(absDuty2, 0);
    } else { // forwards
        motorShield.motorDWrite(absDuty4, 0);
        //motorShield.motorDWrite(absDuty2, 1);
    }             
    prev_current_des4 = current_des4; 

    
    
}

int main (void)
{
    
    // Object for 7th order Cartesian foot trajectory
    BezierCurve rDesFoot_bez_1(2,BEZIER_ORDER_FOOT);
    BezierCurve rDesFoot_bez_2(2,BEZIER_ORDER_FOOT);
    
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

            K_xx                        = input_params[5];    // Foot stiffness N/m
            K_yy                        = input_params[6];    // Foot stiffness N/m
            K_xy                        = input_params[7];    // Foot stiffness N/m
            D_xx                        = input_params[8];    // Foot damping N/(m/s)
            D_yy                        = input_params[9];    // Foot damping N/(m/s)
            D_xy                        = input_params[10];   // Foot damping N/(m/s)
            duty_max                    = input_params[11];   // Maximum duty factor
          
            // Get foot trajectory points
            float foot_pts_1[2*(BEZIER_ORDER_FOOT+1)];
            for(int i = 0; i<2*(BEZIER_ORDER_FOOT+1);i++) {
              foot_pts_1[i] = input_params[12+i];    
            }
            rDesFoot_bez_1.setPoints(foot_pts_1);

            float foot_pts_2[2*(BEZIER_ORDER_FOOT+1)];
            for(int i = 0; i<(BEZIER_ORDER_FOOT+1);i++) {
              foot_pts_2[i] = input_params[13+BEZIER_ORDER_FOOT+i];    
            }
            for(int i = BEZIER_ORDER_FOOT+1; i<2*(BEZIER_ORDER_FOOT+1);i++) {
              foot_pts_2[i] = input_params[11-BEZIER_ORDER_FOOT+i];    
            }
            rDesFoot_bez_2.setPoints(foot_pts_2);
            
            // Attach current loop interrupt
            currentLoop.attach_us(CurrentLoop,current_control_period_us);
                        
            // Setup experiment
            t.reset();
            t.start();
            encoderA.reset();
            encoderB.reset();
            encoderC.reset();
            encoderD.reset();

            motorShield.motorAWrite(0, 0); //turn motor A off
            motorShield.motorBWrite(0, 0); //turn motor B off
            motorShield.motorCWrite(0, 0);
            motorShield.motorDWrite(0, 0);
                         
            // Run experiment
            while( t.read() < start_period + traj_period + end_period) { 
                 
                // Read encoders to get motor states
                angle1 = encoderA.getPulses() *PULSE_TO_RAD + angle1_init;       
                velocity1 = encoderA.getVelocity() * PULSE_TO_RAD;
                 
                angle2 = encoderB.getPulses() * PULSE_TO_RAD + angle2_init;       
                velocity2 = encoderB.getVelocity() * PULSE_TO_RAD;  

                angle3 = encoderC.getPulses() *PULSE_TO_RAD + angle1_init;       
                velocity3 = encoderC.getVelocity() * PULSE_TO_RAD;
                 
                angle3 = encoderD.getPulses() * PULSE_TO_RAD + angle2_init;       
                velocity4 = encoderD.getVelocity() * PULSE_TO_RAD;           
                
                const float th1 = angle1;
                const float th2 = angle2;
                const float dth1= velocity1;
                const float dth2= velocity2;

                const float th3 = angle3;
                const float th4 = angle4;
                const float dth3= velocity3;
                const float dth4= velocity4;
 
                // Calculate the Jacobian
                float Jx_th1 = l_OB * cos(th1) + l_AC * cos(th1 + th2) + l_DE * cos(th1);
                float Jy_th1 = l_OB * sin(th1) + l_AC * sin(th1 + th2) + l_DE * sin(th1);
                float Jx_th2 = l_AC * cos(th1 + th2);
                float Jy_th2 = l_AC * sin(th1 + th2);

                float Jx_th3 = l_OB * cos(th3) + l_AC * cos(th3 + th4) + l_DE * cos(th3);
                float Jy_th3 = l_OB * sin(th3) + l_AC * sin(th3 + th4) + l_DE * sin(th3);
                float Jx_th4 = l_AC * cos(th3 + th4);
                float Jy_th4 = l_AC * sin(th3 + th4);
                                
                // Calculate the forward kinematics (position and velocity)
                float xFoot_1 = l_OB * sin(th1) + l_AC * sin(th1 + th2) + l_DE * sin(th1);
                float yFoot_1 = -l_OB * cos(th1) - l_AC * cos(th1 + th2) - l_DE * cos(th1);
                float dxFoot_1 = l_OB * cos(th1) * dth1 + l_AC * cos(th1 + th2) * (dth1 + dth2) + l_DE * cos(th1) * dth1; 
                float dyFoot_1 = l_OB * sin(th1) * dth1 + l_AC * sin(th1 + th2) * (dth1 + dth2) + l_DE * sin(th1) * dth1;   

                float xFoot_2 = l_OB * sin(th3) + l_AC * sin(th3 + th4) + l_DE * sin(th3);
                float yFoot_2 = -l_OB * cos(th3) - l_AC * cos(th3 + th4) - l_DE * cos(th3);
                float dxFoot_2 = l_OB * cos(th3) * dth3 + l_AC * cos(th3 + th4) * (dth3 + dth4) + l_DE * cos(th3) * dth3; 
                float dyFoot_2 = l_OB * sin(th3) * dth3 + l_AC * sin(th3 + th4) * (dth3 + dth4) + l_DE * sin(th3) * dth3;          

                // Set gains based on buffer and traj times, then calculate desired x,y from Bezier trajectory at current time if necessary
                float teff  = 0;
                float vMult = 0;
                if( t < start_period) {
                    if (K_xx > 0 || K_yy > 0) {
                        K_xx = 100; 
                        K_yy = 100; 
                        D_xx = 5;  
                        D_yy = 5;  
                        K_xy = 0;
                        D_xy = 0;
                    }
                    teff = 0;
                }
                else if (t < start_period + traj_period)
                {
                    K_xx = input_params[5];  // Foot stiffness N/m
                    K_yy = input_params[6];  // Foot stiffness N/m
                    K_xy = input_params[7];  // Foot stiffness N/m
                    D_xx = input_params[8];  // Foot damping N/(m/s)
                    D_yy = input_params[9];  // Foot damping N/(m/s)
                    D_xy = input_params[10]; // Foot damping N/(m/s)
                    teff = (t-start_period);
                    vMult = 1;
                }
                else
                {
                    teff = traj_period;
                    vMult = 0;
                }
                
                // Get desired foot positions and velocities
                float rDesFoot_1[2] , vDesFoot_1[2];
                rDesFoot_bez_1.evaluate(teff/traj_period,rDesFoot_1);
                rDesFoot_bez_1.evaluateDerivative(teff/traj_period,vDesFoot_1);
                vDesFoot_1[0]/=traj_period;
                vDesFoot_1[1]/=traj_period;
                vDesFoot_1[0]*=vMult;
                vDesFoot_1[1]*=vMult;

                float rDesFoot_2[2] , vDesFoot_2[2];
                rDesFoot_bez_2.evaluate(teff/traj_period,rDesFoot_2);
                rDesFoot_bez_2.evaluateDerivative(teff/traj_period,vDesFoot_2);
                vDesFoot_2[0]/=traj_period;
                vDesFoot_2[1]/=traj_period;
                vDesFoot_2[0]*=vMult;
                vDesFoot_2[1]*=vMult;
                
                // Calculate the inverse kinematics (joint positions and velocities) for desired joint angles              
                float xFoot_inv_1 = -rDesFoot_1[0];
                float yFoot_inv_1 = rDesFoot_1[1];                
                float l_OE_1 = sqrt( (pow(xFoot_inv_1,2) + pow(yFoot_inv_1,2)) );
                float alpha_1 = abs(acos( (pow(l_OE_1,2) - pow(l_AC,2) - pow((l_OB+l_DE),2))/(-2.0f*l_AC*(l_OB+l_DE)) ));
                float th2_des_1 = -(3.14159f - alpha_1); 
                float th1_des_1 = -((3.14159f/2.0f) + atan2(yFoot_inv_1,xFoot_inv_1) - abs(asin( (l_AC/l_OE_1)*sin(alpha_1) )));
                
                float dd_1 = (Jx_th1*Jy_th2 - Jx_th2*Jy_th1);
                float dth1_des_1 = (1.0f/dd_1) * (  Jy_th2*vDesFoot_1[0] - Jx_th2*vDesFoot_1[1] );
                float dth2_des_1 = (1.0f/dd_1) * ( -Jy_th1*vDesFoot_1[0] + Jx_th1*vDesFoot_1[1] );

                float xFoot_inv_2 = -rDesFoot_2[0];
                float yFoot_inv_2 = rDesFoot_2[1];                
                float l_OE_2 = sqrt( (pow(xFoot_inv_2,2) + pow(yFoot_inv_2,2)) );
                float alpha_2 = abs(acos( (pow(l_OE_2,2) - pow(l_AC,2) - pow((l_OB+l_DE),2))/(-2.0f*l_AC*(l_OB+l_DE)) ));
                float th2_des_2 = -(3.14159f - alpha_2); 
                float th1_des_2 = -((3.14159f/2.0f) + atan2(yFoot_inv_2,xFoot_inv_2) - abs(asin( (l_AC/l_OE_1)*sin(alpha_2) )));
                
                float dd_2 = (Jx_th1*Jy_th2 - Jx_th2*Jy_th1);
                float dth1_des_2 = (1.0f/dd_2) * (  Jy_th2*vDesFoot_2[0] - Jx_th2*vDesFoot_2[1] );
                float dth2_des_2 = (1.0f/dd_2) * ( -Jy_th1*vDesFoot_2[0] + Jx_th1*vDesFoot_2[1] );
        
                // Calculate error variables
                float e_x_1 = rDesFoot_1[0] - xFoot_1;
                float e_y_1 = rDesFoot_1[1] - yFoot_1;
                float de_x_1 = vDesFoot_1[0] - dxFoot_1;
                float de_y_1 = vDesFoot_1[1] - dyFoot_1;

                float e_x_2 = rDesFoot_2[0] - xFoot_2;
                float e_y_2 = rDesFoot_2[1] - yFoot_2;
                float de_x_2 = vDesFoot_2[0] - dxFoot_2;
                float de_y_2 = vDesFoot_2[1] - dyFoot_2;
        
                // Calculate virtual force on foot
                float fx_1 = K_xx*e_x_1 + K_xy*e_y_1 + D_xx*de_x_1 + D_xy*de_y_1;
                float fy_1 = K_xy*e_x_1 + K_yy*e_y_1 + D_xy*de_x_1 + D_yy*de_y_1;
                //float fx = (K_xx*e_x-K_xx*e_y+2*K_xy*e_x+K_yy*e_x+K_yy*e_y)/2;
                //float fy = (-K_xx*e_x+K_xx*e_y-2*K_xy*e_y+K_yy*e_x+K_yy*e_y)/2;
                
                float fx_2 = K_xx*e_x_2 + K_xy*e_y_2 + D_xx*de_x_2 + D_xy*de_y_2;
                float fy_2 = K_xy*e_x_2 + K_yy*e_y_2 + D_xy*de_x_2 + D_yy*de_y_2;
                
                // Calculate mass matrix elements
                //float M11 = I1 + I2 + I3 + I4 + m1*0.032*0.032 + m2*(l_OB*l_OB+l_B_m2*l_B_m2+2*l_OB*l_B_m2*cos(th2)) + m3*(l_OA*l_OA+l_A_m3*l_A_m3+2*l_OA*l_A_m3*cos(th2)) + m4*((l_OA+l_C_m4)*(l_OA+l_C_m4)+l_AC*l_AC+2*(l_OA+l_C_m4)*l_AC*cos(th2));
                //float M12 = I2 + I3 + m2*(l_B_m2*l_B_m2+l_OB*l_B_m2*cos(th2)) + m3*(l_A_m3*l_A_m3+l_OA*l_A_m3*cos(th2)) + m4*(l_AC*l_AC+(l_OA+l_C_m4)*l_AC*cos(th2));
                //float M22 = I2 + I3 + m2*l_B_m2*l_B_m2 + m3*l_A_m3*l_A_m3 + m4*l_AC*l_AC;
                float t2_1 = cos(th2); float t2_2 = cos(th4);
                float t3 = pow(N,2);
                float t4 = pow(l_AC,2);
                float t5 = pow(l_A_m3,2);
                float t6 = pow(l_B_m2,2);
                float t7 = pow(l_OA,2);
                float t8 = Ir*N;
                float t9 = Ir*t3;
                float t10 = m4*t4;
                float t11 = m3*t5;
                float t12 = m2*t6;
                float t13_1 = l_AC*l_C_m4*m4*t2_1;
                float t14_1 = l_AC*l_OA*m4*t2_1;
                float t15_1 = l_A_m3*l_OA*m3*t2_1;
                float t16_1 = l_B_m2*l_OB*m2*t2_1;
                float t17_1 = I2+I3+t8+t10+t11+t12+t13_1+t14_1+t15_1+t16_1;
                float t13_2 = l_AC*l_C_m4*m4*t2_2;
                float t14_2 = l_AC*l_OA*m4*t2_2;
                float t15_2 = l_A_m3*l_OA*m3*t2_2;
                float t16_2 = l_B_m2*l_OB*m2*t2_2;
                float t17_2 = I2+I3+t8+t10+t11+t12+t13_2+t14_2+t15_2+t16_2;
                //A = reshape([I1+I2+I3+I4+Ir+t9+t10+t11+t12+t13.*2.0+t14.*2.0+t15.*2.0+t16.*2.0+m3.*t7+m4.*t7+l_C_m4.^2.*m4+l_OB.^2.*m2+l_O_m1.^2.*m1+l_C_m4.*l_OA.*m4.*2.0,t17,t17,I2+I3+t9+t10+t11+t12],[2,2]);
                float M11_1 = I1+I2+I3+I4+Ir+t9+t10+t11+t12+t13_1*2.0+t14_1*2.0+t15_1*2.0+t16_1*2.0+m3*t7+m4*t7+pow(l_C_m4,2)*m4+pow(l_OB,2)*m2+pow(l_O_m1,2)*m1+l_C_m4*l_OA*m4*2.0;
                float M12_1 = t17_1;
                float M11_2 = I1+I2+I3+I4+Ir+t9+t10+t11+t12+t13_2*2.0+t14_2*2.0+t15_2*2.0+t16_2*2.0+m3*t7+m4*t7+pow(l_C_m4,2)*m4+pow(l_OB,2)*m2+pow(l_O_m1,2)*m1+l_C_m4*l_OA*m4*2.0;
                float M12_2 = t17_2;
                float M22 = I2+I3+t9+t10+t11+t12;
                
                
                // Populate mass matrix
                MassMatrix_1.Clear();
                MassMatrix_1 << M11_1 << M12_1
                           << M12_1 << M22;

                MassMatrix_2.Clear();
                MassMatrix_2 << M11_2 << M12_2
                           << M12_2 << M22;
                
                // Populate Jacobian matrix
                Jacobian_1.Clear();
                Jacobian_1 << Jx_th1 << Jx_th2
                         << Jy_th1 << Jy_th2;

                Jacobian_2.Clear();
                Jacobian_2 << Jx_th3 << Jx_th4
                         << Jy_th3 << Jy_th4;
                
                // Once you have copied the elements of the mass matrix, uncomment the following section
                
                // Calculate Lambda matrix
                JacobianT_1 = MatrixMath::Transpose(Jacobian_1);
                InverseMassMatrix_1 = MatrixMath::Inv(MassMatrix_1);
                temp_product_1 = Jacobian_1*InverseMassMatrix_1*JacobianT_1;
                Lambda_1 = MatrixMath::Inv(temp_product_1); 

                JacobianT_2 = MatrixMath::Transpose(Jacobian_2);
                InverseMassMatrix_2 = MatrixMath::Inv(MassMatrix_2);
                temp_product_2 = Jacobian_2*InverseMassMatrix_1*JacobianT_2;
                Lambda_2 = MatrixMath::Inv(temp_product_2); 
                
                // Pull elements of Lambda matrix
                float L11_1 = Lambda_1.getNumber(1,1);
                float L12_1 = Lambda_1.getNumber(1,2);
                float L21_1 = Lambda_1.getNumber(2,1);
                float L22_1 = Lambda_1.getNumber(2,2);   

                float L11_2 = Lambda_2.getNumber(1,1);
                float L12_2 = Lambda_2.getNumber(1,2);
                float L21_2 = Lambda_2.getNumber(2,1);
                float L22_2 = Lambda_2.getNumber(2,2);         
                                     
                // Set desired currents             
                current_des1 = ((Jx_th1*L11_1+Jy_th1*L12_1)*fx_1 + (Jx_th1*L12_1+Jy_th1*L22_1)*fy_1) / k_t;          
                current_des2 = ((Jx_th2*L11_1+Jy_th2*L12_1)*fx_1 + (Jx_th2*L12_1+Jy_th2*L22_1)*fy_1) / k_t;   
                current_des3 = ((Jx_th3*L11_2+Jy_th3*L12_2)*fx_2 + (Jx_th3*L12_2+Jy_th3*L22_2)*fy_2) / k_t;          
                current_des4 = ((Jx_th4*L11_2+Jy_th4*L12_2)*fx_2 + (Jx_th4*L12_2+Jy_th4*L22_2)*fy_2) / k_t;   



                // Form output to send to MATLAB     
                float output_data[NUM_OUTPUTS];
                // current time
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
                // foot state
                output_data[11] = xFoot_1;
                output_data[12] = yFoot_1;
                output_data[13] = dxFoot_1;
                output_data[14] = dyFoot_1;
                output_data[15] = rDesFoot_1[0];
                output_data[16] = rDesFoot_1[1];
                output_data[17] = vDesFoot_1[0];
                output_data[18] = vDesFoot_1[1];
                
                // Send data to MATLAB
                server.sendData(output_data,NUM_OUTPUTS);

                wait_us(impedance_control_period_us);   
            }
            
            // Cleanup after experiment
            server.setExperimentComplete();
            currentLoop.detach();
            motorShield.motorAWrite(0, 0); //turn motor A off
            motorShield.motorBWrite(0, 0); //turn motor B off
            motorShield.motorCWrite(0, 0);
            motorShield.motorDWrite(0, 0);
        
        } // end if
        
    } // end while
    
} // end main

