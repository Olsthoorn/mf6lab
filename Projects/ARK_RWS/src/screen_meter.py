import time
from pynput import mouse

class SchermMeter:
    def __init__(self):
        self.stap = 0
        self.y_datapunt = 0
        self.y_as_min = 0
        self.y_as_max = 0
        
        # We vragen eenmalig naar de werkelijke waarden op de Y-as
        print("--- CONFIGURATIE VAN DE Y-AS ---")
        try:
            self.waarde_as_min = float(input("Wat is de onderste waarde op de Y-as? (bijv. 0): ") or 0)
            self.waarde_as_max = float(input("Wat is de bovenste waarde op de Y-as? (bijv. 1): ") or 1)
        except ValueError:
            print("Ongeldige invoer, we gebruiken 0 en 1.")
            self.waarde_as_min = 0.0
            self.waarde_as_max = 1.0
            
        self.schaal_bereik = self.waarde_as_max - self.waarde_as_min
        
        print("\n=== PROGRAMMA GESTART ===")
        print("Ga naar uw rapport en voer de klikken uit.")
        print("STAP 1: Klik op het DATAPUNT op de grafieklijn.")

    def on_click(self, x, y, button, pressed):
        # We luisteren alleen naar het NEERDRUKKEN van de linkermuisknop
        if button == mouse.Button.left and pressed:
            # Op Windows/Mac/Linux is de Y-as (0) altijd de bovenkant van het hoofdscherm.
            # Als uw tweede scherm boven of onder uw hoofdscherm staat, kunnen Y-waarden 
            # negatief of heel groot zijn, maar de onderlinge pixelafstanden blijven exact kloppen.
            
            if self.stap == 0:
                self.y_datapunt = y
                print(f"-> Klik 1 (Datapunt) opgevangen op schermpositie: X={x}, Y={y}")
                print("STAP 2: Klik op de onderste referentiewaarde van de Y-as.")
                self.stap += 1
                
            elif self.stap == 1:
                self.y_as_min = y
                print(f"-> Klik 2 (Y-as minimum) opgevangen op schermpositie: X={x}, Y={y}")
                print("STAP 3: Klik op de bovenste referentiewaarde van de Y-as.")
                self.stap += 1
                
            elif self.stap == 2:
                self.y_as_max = y
                print(f"-> Klik 3 (Y-as maximum) opgevangen op schermpositie: X={x}, Y={y}")
                
                # Berekening van de onderlinge pixelafstanden
                # Omdat de beeldscherm-Y-as omlaag loopt, is het 'Y-minimum' (bijv. 0) 
                # op het scherm juist een HOGER pixelgetal dan het 'Y-maximum' (bijv. 1).
                pixel_afstand_as = self.y_as_min - self.y_as_max
                pixel_afstand_datapunt = self.y_as_min - self.y_datapunt
                
                if pixel_afstand_as == 0:
                    print("\nFout: Klik 2 en Klik 3 liggen op exact dezelfde hoogte! Berekening afgebroken.")
                else:
                    # Lineaire interpolatie inclusief de ingestelde schaal
                    factor = pixel_afstand_datapunt / pixel_afstand_as
                    werkelijke_waarde = self.waarde_as_min + (factor * self.schaal_bereik)
                    
                    print("\n================ BEREKENING ================")
                    print(f"Pixelhoogte totale Y-as op uw scherm : {abs(pixel_afstand_as)} pixels")
                    print(f"Pixelhoogte datapunt vanaf Y-minimum : {pixel_afstand_datapunt} pixels")
                    print(f"Berekende fysieke waarde             : {werkelijke_waarde:.4f}")
                    print("============================================\n")
                    
                # Korte pauze om te voorkomen dat de muis-omhoog actie per ongeluk als nieuwe klik telt
                time.sleep(0.2)
                
                print("STAP 1: Klik opnieuw op een DATAPUNT voor een volgende grafiek.")
                self.stap = 0

# Start de muis-luisteraar
if __name__ == "__main__":
    meter = SchermMeter()
    
    # Dit start een achtergrond-thread die continu naar het hele scherm luistert
    with mouse.Listener(on_click=meter.on_click) as listener:
        listener.join()
