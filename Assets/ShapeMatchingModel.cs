using System.Collections.Generic;
using Unity.Mathematics;
using UnityEngine;
using UnityTools.Common;
using UnityTools.Debuging.EditorTool;
using UnityTools.Math;

namespace UnityClimbingPlant
{
    public class ShapeMatchingModel : MonoBehaviour
    {
        public class Particle : IVertex
        {
            public float3 pos = 0;
            public float3 predictPos = 0;
            public float3 restPos = 0;

            public quaternion rotation = quaternion.identity;
            public quaternion predictRotation = quaternion.identity;

            public float3 w = 0;//angular velocity
            public float3 velocity = 0;

            public float3 ellipsoid = 1;// ellipsoid size
            public float density = 1;
            public object Clone()
            {
                return new Particle() 
                { 
                    pos = pos, predictPos = predictPos, restPos = restPos, 
                    rotation = rotation, predictRotation = predictRotation
                };
            }

            public bool Equals(IVertex other)
            {
                return base.Equals(other);
            }

            public float3x3 Ai
            {
                get
                {
                    var A = new float3x3(this.a * this.a, 0, 0,
                                         0, this.b * this.b, 0,
                                         0, 0, this.c * this.c);

                    var R = new float3x3(this.predictRotation);
                    return (1.0f / 5) * this.Mass * math.mul(A, R);
                }
            }

            public float a => this.ellipsoid.x;
            public float b => this.ellipsoid.y;
            public float c => this.ellipsoid.z;
            public float Mass => this.Volume * this.density;

            public float Volume => 4f / 3f * math.PI * this.a * this.b * this.c;

        }
        public class Edge : IEdge
        {
            public bool IsDirectional => false;

            public IVertex Vertex { get; set; }
            public IVertex OtherVertex { get; set; }

            public object Clone()
            {
                return new Edge() { Vertex = this.Vertex.Clone() as IVertex, OtherVertex = this.OtherVertex.Clone() as IVertex };
            }

            public bool Equals(IEdge other)
            {
                return base.Equals(other);
            }
        }

        public class ParticleGraphFactory : IGraphFactory
        {
            public IEdge CreateEdge(IVertex v1, IVertex v2, bool isDirectional = false)
            {
                return new Edge() { Vertex = v1, OtherVertex = v2 };
            }

            public INewGraph CreateGraph()
            {
                return new ParticleGraph();
            }

            public IVertex CreateVertex()
            {
                return new Particle();
            }
        }

        public class ParticleGraph : NewGraph<Particle, Edge, ParticleGraphFactory>
        {
        }


        [SerializeField] protected Mesh testMesh;
        [SerializeField] protected Transform anchor;
        [SerializeField] protected Bounds domain = new Bounds(new float3(0, 20, 0), new float3(40, 40, 20));

        [SerializeField] protected float stiffness = 1;
        [SerializeField] protected float dt = 0.05f;
        [SerializeField] protected float3 gravity = new float3(0, -9.8f, 0);

        protected ParticleGraph g = new ParticleGraph();
        protected void Start()
        {
            this.AddMesh(this.testMesh);
        }

        protected void Update()
        {
            this.Prediction();
            this.PositionCorrection();
            this.VelocityUpdate();
        }


        protected void AddMesh(Mesh m)
        {
            var transform = this.gameObject.transform;
            var added = new Dictionary<Vector3, Particle>();
            for (var t = 0; t < m.triangles.Length; t += 3)
            {
                var v1 = m.vertices[m.triangles[t]];
                var v2 = m.vertices[m.triangles[t + 1]];
                var v3 = m.vertices[m.triangles[t + 2]];

                Particle p1;
                Particle p2;
                Particle p3;
                if (!added.TryGetValue(v1, out p1))
                {
                    p1 = this.g.AddVertex(this.g.Factory.CreateVertex()) as Particle;
                    added.Add(v1, p1);
                }
                if (!added.TryGetValue(v2, out p2))
                {
                    p2 = this.g.AddVertex(this.g.Factory.CreateVertex()) as Particle;
                    added.Add(v2, p2);
                }
                if (!added.TryGetValue(v3, out p3))
                {
                    p3 = this.g.AddVertex(this.g.Factory.CreateVertex()) as Particle;
                    added.Add(v3, p3);
                }


                p1.pos = p1.restPos = transform.TransformPoint(v1);
                p2.pos = p2.restPos = transform.TransformPoint(v2);
                p3.pos = p3.restPos = transform.TransformPoint(v3);

                this.g.AddEdge(p1, p2);
                this.g.AddEdge(p2, p3);
                this.g.AddEdge(p3, p1);
            }
        }

        protected void Prediction()
        {
            var count = 0;
            foreach (var p in this.g.Vertices)
            {
                p.predictPos = p.pos + p.velocity * dt + 0.5f * (gravity) * dt * dt;
                p.predictPos = math.clamp(p.predictPos, this.domain.min, this.domain.max);
                if (Input.GetKey(KeyCode.A) && count++ < 1) p.predictPos = this.anchor.position;
                
                var w = p.w;
                var wl = math.length(w);

                var predictQ = p.predictRotation;
                if (wl > 0)
                {
                    var wn = math.normalize(w);
                    var wq = wn * math.sin(0.5f * dt * wl);
                    var wdt = new quaternion(wq.x, wq.y, wq.z, math.cos(0.5f * dt * wl));
                    predictQ = math.mul(wdt, p.predictRotation);
                }
                p.predictRotation = predictQ;
            }
        }

        protected void PositionCorrection()
        {
            foreach (var pc in this.g.Vertices)
            {
                this.ApplyDynamicConstraint(pc);
            }
        }

        protected void ApplyDynamicConstraint(Particle pc)
        {
            var plist = new List<Particle>();
            plist.Add(pc);
            foreach (Particle p in this.g.GetNeighborVertices(pc)) plist.Add(p);

            var A = new float3x3(0);//1,0,0,0,1,0,0,0,1);
            foreach (var p in plist)
            {
                var cx = UnityTools.Math.Operation.OuterProduct(p.predictPos, p.restPos);
                var m = p.Mass;
                A += p.Ai + m * cx;
            }
            var C0 = this.C0(plist);
            var Ct = this.Ct(plist);
            var MassSum = this.MassSum(plist);
            var cc = UnityTools.Math.Operation.OuterProduct(Ct, C0);
            A -= MassSum * cc;

            var U = new float3x3(0);
            var d = new float3(0);
            var V = new float3x3(0);

            SVD.GetSVD3D(A, out U, out d, out V);

            var R = math.mul(U, math.transpose(V));


            foreach (var p in plist)
            {
                var goal = math.mul(R, p.restPos - C0) + Ct;
                var delta = goal - p.predictPos;
                
                p.predictPos += delta * this.stiffness;
            }

            pc.predictRotation = new quaternion(R);
        }

        protected void VelocityUpdate()
        {
            foreach (var p in this.g.Vertices)
            {
                p.velocity = (p.predictPos - p.pos) / dt;
                p.pos = p.predictPos;


                var pq = math.mul(p.predictRotation, math.inverse(p.rotation));
                if (pq.value.w < 0) pq.value = -pq.value;

                var angle = 0f;
                var axis = new float3(1, 0, 0);

                this.GetAngleAxis(pq, out angle, out axis);

                if (math.abs(angle) < 0.01f) angle = 0;
                // angle = math.clamp(angle, -1f,1f);
                // p.qt = quaternion.AxisAngle(axis, angle);

                p.w = axis * angle / dt;
                p.rotation = p.predictRotation;
            }
        }

        protected void OnDrawGizmos()
        {
            foreach(var e in this.g.Edges)
            {
                var p1 = e.Vertex as Particle;
                var p2 = e.OtherVertex as Particle;
                Gizmos.DrawLine(p1.pos, p2.pos);
            }

            using(new GizmosScope(Color.green, Matrix4x4.identity))
            {
                Gizmos.DrawWireCube(this.domain.center, this.domain.size);
            }
        }
        public float MassSum(List<Particle> list)
        {
            var ret = 0f;
            foreach (Particle p in list)
            {
                ret += p.Mass;
            }
            return ret;
        }
        public float3 Ct(List<Particle> list)
        {
            var ret = new float3(0);
            var massSum = 0f;
            foreach (Particle p in list)
            {
                ret += p.Mass * p.predictPos;
                massSum += p.Mass;
            }
            return ret / massSum;

        }
        public float3 C(List<Particle> list)
        {
            var ret = new float3(0);
            var massSum = 0f;
            foreach (Particle p in list)
            {
                ret += p.Mass * p.pos;
                massSum += p.Mass;
            }
            return ret / massSum;

        }

        public float3 C0(List<Particle> list)
        {
            var ret = new float3(0);
            var massSum = 0f;
            foreach (Particle p in list)
            {
                ret += p.Mass * p.restPos;
                massSum += p.Mass;
            }
            return ret / massSum;
        }
        protected void GetAngleAxis(quaternion q, out float angle, out float3 axis)
        {
            angle = 0;
            axis = new float3(0, 0, 0);

            var qv = math.normalize(q).value;

            angle = 2 * math.acos(qv.w);
            float s = math.sqrt(1 - qv.w * qv.w);
            if (qv.w == 1)
            {
                axis = qv.xyz;
            }
            else
            {
                axis = qv.xyz / s;
            }
        }
    }
}