#include <SPI.h>
#include <WiFiNINA.h>

// Replace these with your network credentials
const char* ssid = "UOB240020_1634";
const char* password = "mko0)OKM";

int motorPin = 2;

WiFiServer server(80); // Server will run on port 80

void setup() {
  pinMode(motorPin, OUTPUT);
  digitalWrite(motorPin, LOW); // Start with the vibration motors off

  // Initialize serial and wait for the port to open:
  Serial.begin(9600);
  // while (!Serial) {
  //  ; // Wait for serial port to connect.
  // }

  // Connect to WiFi network
  Serial.print("Connecting to ");
  Serial.println(ssid);
  
  WiFi.begin(ssid, password);
  
  while (WiFi.status() != WL_CONNECTED) {
    delay(500);
    Serial.print(".");
  }
  
  Serial.println("");
  Serial.println("WiFi connected.");
  
  // Start the server
  server.begin();
  Serial.println("Server started");

  // Print the IP address
  Serial.print("Use this URL to connect: ");
  Serial.print("http://");
  Serial.print(WiFi.localIP());
  Serial.println("/");
}

void loop() {
  // UNTESTED!!!
  // // Attempt a reconnection if there's a loss of connectivity
  // if (WiFi.status() != WL_CONNECTED) {
  //   Serial.println("Disconnected from WiFi. Attempting reconnection...");
  //   WiFi.begin(ssid, password);
  //   while (WiFi.status() != WL_CONNECTED) {
  //     delay(500);
  //     Serial.print(".");
  //   }
  //   Serial.println("Reconnected.");
  // }

  // Check if a client has connected
  WiFiClient client = server.available();
  if (client) {
    String currentLine = "";
    while (client.connected()) {
      if (client.available()) {
        char c = client.read();
        Serial.write(c);
        if (c == '\n') {
          if (currentLine.length() == 0) {
            // Send a standard HTTP response header
            client.println("HTTP/1.1 200 OK");
            client.println("Content-type:text/html");
            client.println("Connection: close");
            client.println();
            client.println("<!DOCTYPE html><html>");
            client.println("<head><meta name=\"viewport\" content=\"width=device-width, initial-scale=1\"></head>");
            client.println("<body><h1>Arduino Vibration Motor Control</h1>");
            client.println("<p><a href=\"/vibMotor=ON\"><button>ON</button></a></p>");
            client.println("<p><a href=\"/vibMotor=OFF\"><button>OFF</button></a></p>");
            client.println("</body></html>");
            break;
          } else {
            if (currentLine.startsWith("GET /vibMotor=ON")) {
              analogWrite(motorPin, 120);
              // digitalWrite(motorPin, HIGH); // Turn the vibration motor on
            }
            if (currentLine.startsWith("GET /vibMotor=OFF")) {
              analogWrite(motorPin, 0);
              // digitalWrite(motorPin, LOW); // Turn the vibration motor off
            }
            currentLine = "";
          }
        } else if (c != '\r') {
          currentLine += c;
        }
      }
    }
    client.stop();
  }
}
