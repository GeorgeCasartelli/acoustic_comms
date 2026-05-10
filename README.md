# Aerial Acoustic Communications Demonstrator
**Final Year BEng Project**

A MATLAB implementation of an OFDM modem for data transmission via acoustic waves. This system is designed to transmit text data or small images through air using 400 subcarriers with pilot based channel correction.

## System Overview
The project implements a digital communication pipeline adapted for acoustic frequencies. Preamble correlation allows for packets to be detected. Pilot symbols allow for equalisation of the acoustic channel, fixing Carrier Frequency Offset (CFO). Synchronisation occurs with known preamble based correlation, and message length is embedded in packet header.

## Possible implementations

Plan is to use system as a method of generating audio QR codes. Phone app will be built as a receiver but same tech will be used, maybe converted to python app but who knows.

## Hardware Setup
- **Transmitter:** Genelec Active Studio Monitor
- **Receiver:** Shure SM57 Dynamic Microphone
- **Interfaces:** Two Focusrite Scarlett interfaces (Asynchronous)

## Prerequisites
- MATLAB
- Communications Toolbox
- Audio Toolbox
- Signal Processing Toolbox

## Execution
Current WIP for OFDM.

## Placeholder Metrics
| Parameter | Value |
| :--- | :--- |
| **Bandwidth** | 9375 kHz |
| **Bit Rate** | 4997 bps |
| **Target SNR** | [Placeholder] dB |
| **Effective Range** | 2+ m |
| **Carrier Freq** | 10000 Hz |
| **Samplerate** | 48000 Hz |
