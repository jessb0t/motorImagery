## Motor power modulations during imagined movements

#### Abstract

Neural activity patterns shift in human motor cortex before and during voluntary movements. Prior work has demonstrated a decrease in relative power in the beta frequency band during sustained movement. Likewise, an increase in relative power in the high gamma frequency band has been found to correlate with the dynamic features of a movement, such as individual finger flexion. These rhythms may indicate the presence of neuronal population dynamics involved in initiating movement.

Despite the widespread use of imagined motor movements for control of brain-computer interfaces, few studies have addressed the neuronal basis of imagined movements. We hypothesized that imagined movement would demonstrate amplitude modulations similar to actual movement across relevant motor areas, albeit with an attenuated response in higher frequency ranges given the absence of actual movement.

After confirming this hypothesis, we conjectured that these power modulations could be relevant features for training a classifier algorithm to identify imagined versus actual movements. Although we believed that high gamma would be the most relevant differentiator, not all imaging techniques have access to this high power band. We therefore investigated whether a decent classifier could be achieved by neglecting high gamma in favor of lower frequency bands.

We reanalyzed electrocortigraphy (ECoG) recordings from Miller et al. (2019) from subjects undergoing monitoring for treatment of medically refractory epilepsy. Subjects performed two tasks involving actual or imagined movement of the hand or tongue. Dataglove or electromyography measurement verified an absence of motor action during imagined movement. Following Unterweger et al. (2020), we defined six functionally relevant frequency bands from 12-150 Hz. Event related synchrony and desynchronization (ERS/ERD) measurements were determined by relative amplitude changes in respective power spectra during movement initiation.

Using Random Forest models, we compared classification performance with and without high gamma. Although models that included high gamma performed better, we achieved 88% ROC with lower frequency bands alone, indicating that lower frequencies substantially contribute to the distinction between actual and imagined motor movements.

Power modulations appear to serve as a measure of movement intent, which can be used as improved control signals for brain-computer interface applications. Knowledge of these dynamics may additionally enhance our knowledge of subcortical-cortical and cortico-cortical communication in motor control to serve as a biomarker in health and disease.


#### Summary
This group project was completed during Neuromatch 2022 (July 11-29). See `project-summary.pdf` for details!

#### Contributors
Brenda Qiu, Jessica Alexander, Juan Pablo Botero, Kurt Lehner, Lavanya M K