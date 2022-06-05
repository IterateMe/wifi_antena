SYSTEM_THREAD(ENABLED);
SYSTEM_MODE(SEMI_AUTOMATIC);

void setup() {
  WiFi.on();
  Serial.begin(9600);
  // Wait for a USB serial connection for up to 30 seconds
  waitFor(Serial.isConnected, 30000);
  if(Serial.isConnected() == false && Particle.connected() == false) {
        Particle.connect();
  }

}


void wifi_scan_callback(WiFiAccessPoint* wap, void* data)
{
    WiFiAccessPoint& ap = *wap;
    Serial.printlnf("ssid=%s security=%d channel=%d rssi=%d", ap.ssid , (int)ap.security, (int)ap.channel, ap.rssi);
    
}

void loop()
{
    int result_count = WiFi.scan(wifi_scan_callback);
    Serial.printlnf("Scan complete found %d APs.", result_count);

}
