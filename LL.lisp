;(ql:quickload 'lispbuilder-sdl)
;(ql:quickload 'lispbuilder-sdl-gfx)

(defparameter *width* 100)
(defparameter *height* 100)
(defparameter *pixsize* 4) ;;Note + TODO: One "pixel" in the sim is equal to 20um, this is an estimate, revise in stage two.
(defparameter *umperpix* 5)

(defparameter *cells* nil)
(defparameter *world* nil)
(defparameter *datafile* "Celldata.db")

(defun random-range (min max)
  (+ min (random (- max min) (make-random-state t))))

(defun new-cell (POS ATP NA AA FA G O2 CO2 DNA)
  (push (list :POS POS :ATP ATP :NA NA 
         :AA AA :FA FA :G G :O2 O2 :CO2 CO2 :DNA DNA) *cells*))

(defun new-world (TEMP PH NAC AAC FAC GC O2C CO2C RAD BPV)
  (setf *world* (list :TEMP TEMP :PH PH :NAC NAC 
                :AAC AAC :FAC FAC :GC GC :O2C O2C :CO2C CO2C :RAD RAD :BPV BPV)))

(defmacro fetch-value (accessor index lst)
  `(getf (nth ,index ,lst) ,accessor))
  
(defun fetch-props (lst)
  (let ((nlst nil))
    (do ((i 0 (+ i 2)))
        ((>= i (length lst)))
      (push (nth i lst) nlst)) (reverse nlst)))  

(defun parse-gene (gene)
  (let ((gene-val 0))
    (dotimes (i (length gene))
      (setf gene-val (+ (* (expt 2 i) (nth i (reverse gene))) gene-val))) gene-val))
      
(defun parse-DNA (DNA-seq)
  (let ((DNA-val nil))
    (dotimes (i (length DNA-seq))
      (push (parse-gene (nth i DNA-seq)) DNA-val)) (reverse DNA-val)))

(defun dump-ci ()
  (dolist (entry *cells*)
    (format t "~%~{~a:~10t~a~%~}~%" entry)))

(defun save-ci (filename)
  (with-open-file (out filename
                   :direction :output
                   :if-exists :supersede)
    (with-standard-io-syntax
      (print *cells* out)
      (print *world* out))))

(defun load-ci (filename)
  (with-open-file (in filename)
    (with-standard-io-syntax
      (setf *cells* (read in))
      (setf *world* (read in)))))
      
(defun rand-mutate (DNA-seq)
  (let ((gene (random-range 0 (length DNA-seq))))
    (let ((base (random-range 0 (length (nth gene DNA-seq)))))
      (cond ((= (nth base (nth gene DNA-seq)) 0) (setf (nth base (nth gene DNA-seq)) 1))
            (t (setf (nth base (nth gene DNA-seq)) 0))) DNA-seq)))
            
(defun mutate-times (DNA-seq times)
  (dotimes (i times)
    (rand-mutate DNA-seq)))
  
(defun close-sim ()
  (format t "Saving simulation state and exiting...~%")
  (save-ci *datafile*)
  (exit))
  
(defun init-sim ()
  (format t "Loading simulation state and starting...~%")
  (when (probe-file *datafile*)
    (load-ci *datafile*))
  (format t "Done."))

;;TODO: Add capability to expel molicules as well.
(defun Cell-Env (celli) ;;Note + TODO: This is scientifically flawed, revise in stage two...
  (let ((prop nil) (cprop nil) (propcont nil) (cpropcont nil) (cellperinc nil) (area (* *width* *height*)))
    (do ((i 0 (1+ i)) (chance (random-range 0 100)))
        ((>= i 6)) ;;Note + TODO: 6 is for the four main macromolicules plus CO2 and O2, make more portable later...
      
      (setf prop (nth i '(:NAC :AAC :FAC :GC :O2C :CO2C))) ;;TODO: Make more portable later...
      (setf cprop (nth i '(:NA :AA :FA :G :O2 :CO2))) ;;TODO: Make more portable later...
      (setf propcont (getf *world* prop))
      (setf cpropcont (fetch-value cprop celli *cells*))
      (setf cellperinc (float (/ propcont (* 5 area)))) ;;Dummy Value: 5 is a placeholder for permeability. <<-- TODO: Revise this
      
      (cond ((<= chance (/ propcont area)) (setf (fetch-value cprop celli *cells*) (+ cpropcont cellperinc)) (setf (getf *world* prop) (- propcont cellperinc)))
            (t (continue))))))
  
(defun Cell-Meta (celli) ;;TODO: Add a way to regulate the ATP synthase
  (cond ((> 1 0) ;;Note + TODO: This is the condition for mitochondrial resperation, revise during stage two...
    (dotimes (i 5) ;;Dummy Value: 5 is a placeholder for the mitochondron count.
     (cond ((and (>= (fetch-value :O2 celli *cells*) 6) (> (fetch-value :G celli *cells*) 0)) ;;Note: Aerobic Resperation C6H12O6 + 6 O2 → 6 CO2 + 6 H2O + 38 ATP. For now, Water is ignored
           (decf (fetch-value :G celli *cells*) 1) (decf (fetch-value :O2 celli *cells*) 6)
           (incf (fetch-value :CO2 celli *cells*) 6) (incf (fetch-value :ATP celli *cells*) 38))
           
           ((and (< (fetch-value :O2 celli *cells*) 6) (> (fetch-value :G celli *cells*) 0)) ;;Note: Anaerobic Resperation C6H12O6 → 2 CO2 + 2 C2H5OH + 2 ATP. For now Ethanol is ignored
           (decf (fetch-value :G celli *cells*) 1) (incf (fetch-value :CO2 celli *cells*) 2)
           (incf (fetch-value :ATP celli *cells*) 2)))))
        ;;TODO: This is the condition for ATP synthase via a proton gradient.   
        (t (format t "If you see this messsage, the laws of science have broken down.~%Now is the time for panic"))))  
  
(defun next-tick ()
  (dotimes (i (length *cells*))
    (Cell-Env i)
    (Cell-Meta i)))

;(defun display-sim (&optional delay)
;  (sdl:with-init ()
;  (sdl:window (* *pixsize* *width*) (* *pixsize* *height*) :title-caption "Lisp Land")
;    (setf (sdl:frame-rate) 60)
;    (sdl:with-events ()
;      (:quit-event () t)
;      (:idle ()
;        (sleep delay)
;        (next-tick)
;        (dotimes (i (length *cells*))
;          (let ((POS (fetch-value :POS i *cells*)))
;                  (sdl-gfx:draw-box (sdl:rectangle :x (* (car POS) *pixsize*) :y (* (cdr POS) *pixsize*)
;                                    :w *pixsize* :h *pixsize*) :color sdl:*white*)))
;        (sdl:update-display)))))

(init-sim)
