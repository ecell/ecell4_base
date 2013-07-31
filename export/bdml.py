from export import Exporter
import xml.etree.ElementTree as ET
from export import particle_spatiocyte_loader

class BdmlExporter(Exporter):
    def export(self):
        component = ET.Element('component')        
        for i in range(181):
            time = ET.SubElement(component, 'time')
            time.text = str(i*.5)
            
            hoge = particle_spatiocyte_loader.load_particles_from_spatiocyte(self.model, index=i)

            measurementList = ET.SubElement(component, 'measurementList')

            for s in hoge.species:
                measurement = ET.SubElement(measurementList, 'measurement')
                targetref = ET.SubElement(measurement, 'targetRef')
                targetref.text = s
                for particle in hoge.list_particles(sid=s):
                    point = ET.SubElement(targetref, 'point')
                    point.text = ','.join(map(str, particle[1].position))

        f = open(self.model.replace(".dat", ".bdml"), 'w')
        f.write(ET.tostring(component))
        f.close()

