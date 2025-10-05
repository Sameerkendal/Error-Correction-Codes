%% parameters

n = 10;  
k = 5;    
tc = floor((n - k)/2);  
prim_poly = 285;  
ms = [23,26,69,244,254];

%% RS ENCODING

m = gf(ms, 8, prim_poly);   
[G_sys,G] = get_systematic_generator_matrix(k, n, 8);

c = m * G_sys;  

%% CHANNEL TRANSMISSION
 Errors = [0,0,2,0,0,3,0,0,0,0];  
 channelErrors = gf(Errors, 8, prim_poly);  
% 
y = c + channelErrors;   

%y = [77	147	248	122	103	19	244	140	217	94	182	220	205	54	165	182	15	147	229	41	123	78	221	163	111	105	186	249	219	87	125	177	36	193	205	50	148	236	49	68	39	26	60	168	177	197	15	37	247	182	109	167	216	26	186	126	101	29	57	119	46	28	222	42	147	95	113	96	65	235	179	33	239	16	108	235	217	20	1	20	40	2	219	101	142	98	128	155	6	77	71	144	70	130	104	188	150	255	251	180	100	107	57	47	195	135	106	236	145	91	216	165	203	106	253	51	60	248	248	58	52	0	36	200	234];
%% DECODING

evalpts = [0:n-1];%[0,1,2,3,4,5,6,7,8,9];
evalpts = gf(evalpts, 8, prim_poly); 
decoded_msg = berlekamp_welch_decoder(evalpts,y,n, k, tc, prim_poly);
disp(decoded_msg.x);


